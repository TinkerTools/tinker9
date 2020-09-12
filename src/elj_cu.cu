#include "add.h"
#include "evdw.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "seq_pair_lj.h"
#include "seq_switch.h"
#include "seq_triangle.h"
#include "spatial2.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
template <class Ver>
__global__
void elj_cu1(count_buffer restrict nebuf, energy_buffer restrict ebuf,
             virial_buffer restrict vbuf, grad_prec* restrict gx,
             grad_prec* restrict gy, grad_prec* restrict gz,
             TINKER_IMAGE_PARAMS, real cut, real off, const real* restrict x,
             const real* restrict y, const real* restrict z, int n,
             const Spatial::SortedAtom* restrict sorted, int nakpl,
             const int* restrict iakpl, int niak, const int* restrict iak,
             const int* restrict lst, int nexclude,
             const int (*restrict exclude)[2],
             const real* restrict exclude_scale,
             const unsigned int* restrict info, int njvdw,
             const real* restrict radmin, const real* restrict epsilon,
             const int* restrict jvdw)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   using ebuf_prec = energy_buffer_traits::type;
   using vbuf_prec = virial_buffer_traits::type;
   ebuf_prec etl;
   vbuf_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
   if CONSTEXPR (do_e) {
      etl = 0;
   }
   if CONSTEXPR (do_v) {
      vtlxx = 0;
      vtlyx = 0;
      vtlzx = 0;
      vtlyy = 0;
      vtlzy = 0;
      vtlzz = 0;
   }
   int ctl;
   if CONSTEXPR (do_a) {
      ctl = 0;
   }


   //* /
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scale = exclude_scale[ii];


      int ijvdw = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      int kjvdw = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real r2 = image2(xr, yr, zr);
      real r = REAL_SQRT(r2);
      real invr = REAL_RECIP(r);
      real rv = radmin[ijvdw * njvdw + kjvdw];
      real eps = epsilon[ijvdw * njvdw + kjvdw];
      real e, de;
      pair_lj_v3<do_g, 0>(r, invr, scale, rv, eps, cut, off, e, de);


      if CONSTEXPR (do_e) {
         etl += cvt_to<ebuf_prec>(e);
         if CONSTEXPR (do_a) {
            if (e != 0) {
               ctl += 1;
            }
         }
      }
      if CONSTEXPR (do_g) {
         real dedx, dedy, dedz;
         de *= invr;
         dedx = de * xr;
         dedy = de * yr;
         dedz = de * zr;
         atomic_add(dedx, gx, i);
         atomic_add(dedy, gy, i);
         atomic_add(dedz, gz, i);
         atomic_add(-dedx, gx, k);
         atomic_add(-dedy, gy, k);
         atomic_add(-dedz, gz, k);
         if CONSTEXPR (do_v) {
            vtlxx += cvt_to<vbuf_prec>(xr * dedx);
            vtlyx += cvt_to<vbuf_prec>(yr * dedx);
            vtlzx += cvt_to<vbuf_prec>(zr * dedx);
            vtlyy += cvt_to<vbuf_prec>(yr * dedy);
            vtlzy += cvt_to<vbuf_prec>(zr * dedy);
            vtlzz += cvt_to<vbuf_prec>(zr * dedz);
         }
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      real fix, fiy, fiz;
      real fkx, fky, fkz;
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }


      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);


      int shiid = ty * WARP_SIZE + ilane;
      int shatomi = min(shiid, n - 1);
      real shxi = sorted[shatomi].x;
      real shyi = sorted[shatomi].y;
      real shzi = sorted[shatomi].z;
      int shi = sorted[shatomi].unsorted;
      int shijvdw = jvdw[shi];


      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      real xk = sorted[atomk].x;
      real yk = sorted[atomk].y;
      real zk = sorted[atomk].z;
      int k = sorted[atomk].unsorted;
      int kjvdw = jvdw[k];


      int bit0 = info[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int srcmask = 1 << srclane;
         int bit = bit0 & srcmask;
         int ijvdw = shijvdw;
         int iid = shiid;
         real xr = shxi - xk;
         real yr = shyi - yk;
         real zr = shzi - zk;


         bool incl = iid < kid and kid < n and bit == 0;
         real r2 = image2(xr, yr, zr);
         real r = REAL_SQRT(r2);
         real invr = REAL_RECIP(r);
         real rv = radmin[ijvdw * njvdw + kjvdw];
         real eps = epsilon[ijvdw * njvdw + kjvdw];
         real e, de;
         pair_lj_v3<do_g, 1>(r, invr, 1, rv, eps, cut, off, e, de);


         if CONSTEXPR (do_e) {
            etl += incl ? cvt_to<ebuf_prec>(e) : 0;
            if CONSTEXPR (do_a) {
               if (incl and e != 0) {
                  ctl += 1;
               }
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de = incl ? de * invr : 0;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            fix += dedx;
            fiy += dedy;
            fiz += dedz;
            fkx -= dedx;
            fky -= dedy;
            fkz -= dedz;
            if CONSTEXPR (do_v) {
               vtlxx += cvt_to<vbuf_prec>(xr * dedx);
               vtlyx += cvt_to<vbuf_prec>(yr * dedx);
               vtlzx += cvt_to<vbuf_prec>(zr * dedx);
               vtlyy += cvt_to<vbuf_prec>(yr * dedy);
               vtlzy += cvt_to<vbuf_prec>(zr * dedy);
               vtlzz += cvt_to<vbuf_prec>(zr * dedz);
            }
         }


         shijvdw = __shfl_sync(ALL_LANES, shijvdw, ilane + 1);
         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
         shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
         shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
         shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
         fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
         fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
         fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, shi);
         atomic_add(fiy, gy, shi);
         atomic_add(fiz, gz, shi);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   } // end loop block pairs
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real fix, fiy, fiz;
      real fkx, fky, fkz;
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }


      int ty = iak[iw];
      int shatomi = ty * WARP_SIZE + ilane;
      real shxi = sorted[shatomi].x;
      real shyi = sorted[shatomi].y;
      real shzi = sorted[shatomi].z;
      int shi = sorted[shatomi].unsorted;
      int shijvdw = jvdw[shi];


      int atomk = lst[iw * WARP_SIZE + ilane];
      real xk = sorted[atomk].x;
      real yk = sorted[atomk].y;
      real zk = sorted[atomk].z;
      int k = sorted[atomk].unsorted;
      int kjvdw = jvdw[k];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int ijvdw = shijvdw;
         real xr = shxi - xk;
         real yr = shyi - yk;
         real zr = shzi - zk;


         bool incl = atomk > 0;
         real r2 = image2(xr, yr, zr);
         real r = REAL_SQRT(r2);
         real invr = REAL_RECIP(r);
         real rv = radmin[ijvdw * njvdw + kjvdw];
         real eps = epsilon[ijvdw * njvdw + kjvdw];
         real e, de;
         pair_lj_v3<do_g, 1>(r, invr, 1, rv, eps, cut, off, e, de);


         if CONSTEXPR (do_e) {
            etl += incl ? cvt_to<ebuf_prec>(e) : 0;
            if CONSTEXPR (do_a) {
               if (incl and e != 0) {
                  ctl += 1;
               }
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de = incl ? de * invr : 0;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            fix += dedx;
            fiy += dedy;
            fiz += dedz;
            fkx -= dedx;
            fky -= dedy;
            fkz -= dedz;
            if CONSTEXPR (do_v) {
               vtlxx += cvt_to<vbuf_prec>(xr * dedx);
               vtlyx += cvt_to<vbuf_prec>(yr * dedx);
               vtlzx += cvt_to<vbuf_prec>(zr * dedx);
               vtlyy += cvt_to<vbuf_prec>(yr * dedy);
               vtlzy += cvt_to<vbuf_prec>(zr * dedy);
               vtlzz += cvt_to<vbuf_prec>(zr * dedz);
            }
         }


         shijvdw = __shfl_sync(ALL_LANES, shijvdw, ilane + 1);
         shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
         shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
         shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
         fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
         fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
         fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, shi);
         atomic_add(fiy, gy, shi);
         atomic_add(fiz, gz, shi);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   } // end loop block-atoms
   // */


   if CONSTEXPR (do_a) {
      atomic_add(ctl, nebuf, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(etl, ebuf, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
   }
}


// special vdw14 interactions
template <class Ver>
__global__
void elj_cu2(count_buffer restrict nebuf, energy_buffer restrict ebuf,
             virial_buffer restrict vbuf, grad_prec* restrict gx,
             grad_prec* restrict gy, grad_prec* restrict gz,
             TINKER_IMAGE_PARAMS, real cut, real off, const real* restrict x,
             const real* restrict y, const real* restrict z, //
             int njvdw, const real* restrict radmin,
             const real* restrict epsilon,
             const int* restrict jvdw, //
             real v4scale, int nvdw14, const int (*restrict vdw14ik)[2],
             const real* restrict radmin4, const real* restrict epsilon4)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   using ebuf_prec = energy_buffer_traits::type;
   using vbuf_prec = virial_buffer_traits::type;
   ebuf_prec etl;
   vbuf_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
   if CONSTEXPR (do_e) {
      etl = 0;
   }
   if CONSTEXPR (do_v) {
      vtlxx = 0;
      vtlyx = 0;
      vtlzx = 0;
      vtlyy = 0;
      vtlzy = 0;
      vtlzz = 0;
   }


   int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   for (int ii = ithread; ii < nvdw14; ii += blockDim.x * gridDim.x) {
      int i = vdw14ik[ii][0];
      int k = vdw14ik[ii][1];


      int jit = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      int jkt = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      int pos = jit * njvdw + jkt;
      real rv = radmin[pos];
      real eps = epsilon[pos];
      real rv4 = radmin4[pos];
      real eps4 = epsilon4[pos];


      real r2 = image2(xr, yr, zr);
      real r = REAL_SQRT(r2);
      real invr = REAL_RECIP(r);


      real e, de, e4, de4;
      pair_lj_v3<do_g, 0>(r, invr, v4scale, rv, eps, cut, off, e, de);
      pair_lj_v3<do_g, 0>(r, invr, v4scale, rv4, eps4, cut, off, e4, de4);
      e = e4 - e;
      if CONSTEXPR (do_g)
         de = de4 - de;


      if CONSTEXPR (do_e) {
         // if CONSTEXPR (do_a) {}
         etl += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_g) {
         de *= invr;
         real dedx = de * xr;
         real dedy = de * yr;
         real dedz = de * zr;
         if CONSTEXPR (do_v) {
            vtlxx += cvt_to<vbuf_prec>(xr * dedx);
            vtlyx += cvt_to<vbuf_prec>(yr * dedx);
            vtlzx += cvt_to<vbuf_prec>(zr * dedx);
            vtlyy += cvt_to<vbuf_prec>(yr * dedy);
            vtlzy += cvt_to<vbuf_prec>(zr * dedy);
            vtlzz += cvt_to<vbuf_prec>(zr * dedz);
         }
         atomic_add(dedx, gx, i);
         atomic_add(dedy, gy, i);
         atomic_add(dedz, gz, i);
         atomic_add(-dedx, gx, k);
         atomic_add(-dedy, gy, k);
         atomic_add(-dedz, gz, k);
      }
   }


   if CONSTEXPR (do_e) {
      atomic_add(etl, ebuf, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
   }
}


template <class Ver>
void elj_cu4()
{
   const auto& st = *cspatial_v2_unit;
   const real cut = switch_cut(switch_vdw);
   const real off = switch_off(switch_vdw);


   int ngrid = get_grid_size(BLOCK_DIM);
   elj_cu1<Ver><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      nev, ev, vir_ev, devx, devy, devz, TINKER_IMAGE_ARGS, cut, off, st.x,
      st.y, st.z, //
      n, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, nvexclude,
      vexclude, vexclude_scale, st.si2.bit0, //
      njvdw, radmin, epsilon, jvdw);
}


template <class Ver>
void elj_cu5()
{
   if (nvdw14 > 0) {
      const auto& st = *cspatial_v2_unit;
      const real cut = switch_cut(switch_vdw);
      const real off = switch_off(switch_vdw);


      launch_k1s(nonblk, nvdw14, elj_cu2<Ver>, //
                 nev, ev, vir_ev, devx, devy, devz, TINKER_IMAGE_ARGS, cut, off,
                 st.x, st.y, st.z,             //
                 njvdw, radmin, epsilon, jvdw, //
                 v4scale, nvdw14, vdw14ik, radmin4, epsilon4);
   }
}


void elj14_cu(int vers)
{
   if (vers == calc::v0)
      elj_cu5<calc::V0>();
   else if (vers == calc::v1)
      elj_cu5<calc::V1>();
   else if (vers == calc::v3)
      elj_cu5<calc::V3>();
   else if (vers == calc::v4)
      elj_cu5<calc::V4>();
   else if (vers == calc::v5)
      elj_cu5<calc::V5>();
   else if (vers == calc::v6)
      elj_cu5<calc::V6>();
}


void elj_cu(int vers)
{
   if (vers == calc::v0)
      elj_cu4<calc::V0>();
   else if (vers == calc::v1)
      elj_cu4<calc::V1>();
   else if (vers == calc::v3)
      elj_cu4<calc::V3>();
   else if (vers == calc::v4)
      elj_cu4<calc::V4>();
   else if (vers == calc::v5)
      elj_cu4<calc::V5>();
   else if (vers == calc::v6)
      elj_cu4<calc::V6>();


   elj14_cu(vers);
}
}
