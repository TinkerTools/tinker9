#include "add.h"
#include "ff/hippo/echgtrn.h"
#include "ff/image.h"
#include "ff/switch.h"
#include "glob/spatial.h"
#include "launch.h"
#include "md/inc.h"
#include "mod/elecpchg.h"
#include "seq/bsplgen.h"
#include "seq/pair_chgtrn.h"
#include "seq/switch.h"
#include "seq/triangle.h"
#include "tool/gpucard.h"
#include <cassert>

namespace tinker {
// ck.py Version 2.0.3
template <class Ver>
__global__
void echgtrn_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nc, energy_buffer restrict ec,
   virial_buffer restrict vc, grad_prec* restrict gx, grad_prec* restrict gy,
   grad_prec* restrict gz, real cut, real off, const unsigned* restrict minfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, real* restrict chgct, real* restrict dmpct,
   real f)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int nctl;
   if CONSTEXPR (do_a) {
      nctl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec ectl;
   if CONSTEXPR (do_e) {
      ectl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz;
   if CONSTEXPR (do_v) {
      vctlxx = 0;
      vctlyx = 0;
      vctlzx = 0;
      vctlyy = 0;
      vctlzy = 0;
      vctlzz = 0;
   }
   real xi;
   real yi;
   real zi;
   real xk;
   real yk;
   real zk;
   real gxi;
   real gyi;
   real gzi;
   real gxk;
   real gyk;
   real gzk;
   real chgi;
   real alphai;
   real chgk;
   real alphak;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii][0];

      xi = x[i];
      yi = y[i];
      zi = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      chgi = chgct[i];
      alphai = dmpct[i];
      chgk = chgct[k];
      alphak = dmpct[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         e_prec e, de;
         pair_chgtrn<do_g>(r, cut, off, scalea, f, alphai, chgi, alphak, chgk, e, de);
         if CONSTEXPR (do_a)
            if (e != 0 and scalea != 0)
               nctl += 1;
         if CONSTEXPR (do_e)
            ectl += floatTo<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(r);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;

            gxi -= dedx;
            gyi -= dedy;
            gzi -= dedz;
            gxk += dedx;
            gyk += dedy;
            gzk += dedz;

            if CONSTEXPR (do_v) {
               vctlxx += floatTo<vbuf_prec>(xr * dedx);
               vctlyx += floatTo<vbuf_prec>(yr * dedx);
               vctlzx += floatTo<vbuf_prec>(zr * dedx);
               vctlyy += floatTo<vbuf_prec>(yr * dedy);
               vctlzy += floatTo<vbuf_prec>(zr * dedy);
               vctlzz += floatTo<vbuf_prec>(zr * dedz);
            }
         }
      } // end if (include)

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }

      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);

      int iid = ty * WARP_SIZE + ilane;
      int atomi = min(iid, n - 1);
      int i = sorted[atomi].unsorted;
      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      int k = sorted[atomk].unsorted;
      xi = sorted[atomi].x;
      yi = sorted[atomi].y;
      zi = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      chgi = chgct[i];
      alphai = dmpct[i];
      chgk = chgct[k];
      alphak = dmpct[k];

      unsigned int minfo0 = minfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (minfo0 & srcmask) == 0;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            e_prec e, de;
            pair_chgtrn<do_g>(r, cut, off, 1, f, alphai, chgi, alphak, chgk, e, de);
            if CONSTEXPR (do_a)
               if (e != 0)
                  nctl += 1;
            if CONSTEXPR (do_e)
               ectl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(r);
               real dedx = de * xr;
               real dedy = de * yr;
               real dedz = de * zr;

               gxi -= dedx;
               gyi -= dedy;
               gzi -= dedz;
               gxk += dedx;
               gyk += dedy;
               gzk += dedz;

               if CONSTEXPR (do_v) {
                  vctlxx += floatTo<vbuf_prec>(xr * dedx);
                  vctlyx += floatTo<vbuf_prec>(yr * dedx);
                  vctlzx += floatTo<vbuf_prec>(zr * dedx);
                  vctlyy += floatTo<vbuf_prec>(yr * dedy);
                  vctlzy += floatTo<vbuf_prec>(zr * dedy);
                  vctlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         chgi = __shfl_sync(ALL_LANES, chgi, ilane + 1);
         alphai = __shfl_sync(ALL_LANES, alphai, ilane + 1);
         if CONSTEXPR (do_g) {
            gxi = __shfl_sync(ALL_LANES, gxi, ilane + 1);
            gyi = __shfl_sync(ALL_LANES, gyi, ilane + 1);
            gzi = __shfl_sync(ALL_LANES, gzi, ilane + 1);
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }

      int ty = iak[iw];
      int atomi = ty * WARP_SIZE + ilane;
      int i = sorted[atomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      xi = sorted[atomi].x;
      yi = sorted[atomi].y;
      zi = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      chgi = chgct[i];
      alphai = dmpct[i];
      chgk = chgct[k];
      alphak = dmpct[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         bool incl = atomk > 0;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            e_prec e, de;
            pair_chgtrn<do_g>(r, cut, off, 1, f, alphai, chgi, alphak, chgk, e, de);
            if CONSTEXPR (do_a)
               if (e != 0)
                  nctl += 1;
            if CONSTEXPR (do_e)
               ectl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(r);
               real dedx = de * xr;
               real dedy = de * yr;
               real dedz = de * zr;

               gxi -= dedx;
               gyi -= dedy;
               gzi -= dedz;
               gxk += dedx;
               gyk += dedy;
               gzk += dedz;

               if CONSTEXPR (do_v) {
                  vctlxx += floatTo<vbuf_prec>(xr * dedx);
                  vctlyx += floatTo<vbuf_prec>(yr * dedx);
                  vctlzx += floatTo<vbuf_prec>(zr * dedx);
                  vctlyy += floatTo<vbuf_prec>(yr * dedy);
                  vctlzy += floatTo<vbuf_prec>(zr * dedy);
                  vctlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)

         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         chgi = __shfl_sync(ALL_LANES, chgi, ilane + 1);
         alphai = __shfl_sync(ALL_LANES, alphai, ilane + 1);
         if CONSTEXPR (do_g) {
            gxi = __shfl_sync(ALL_LANES, gxi, ilane + 1);
            gyi = __shfl_sync(ALL_LANES, gyi, ilane + 1);
            gzi = __shfl_sync(ALL_LANES, gzi, ilane + 1);
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(nctl, nc, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(ectl, ec, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz, vc, ithread);
   }
}

template <class Ver>
void echgtrn_cu2()
{
   const auto& st = *mspatial_v2_unit;
   real cut = switch_cut(switch_chgtrn);
   real off = switch_off(switch_chgtrn);
   real f = electric / dielec;

   assert(ctrntyp == chgtrn_t::SEPARATE);
   int ngrid = gpuGridSize(BLOCK_DIM);
   echgtrn_cu1<Ver><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nct, ect, vir_ect,
      dectx, decty, dectz, cut, off, st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale, st.x,
      st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, chgct, dmpct, f);
}

void echgtrn_cu(int vers)
{
   if (vers == calc::v0)
      echgtrn_cu2<calc::V0>();
   else if (vers == calc::v1)
      echgtrn_cu2<calc::V1>();
   else if (vers == calc::v3)
      echgtrn_cu2<calc::V3>();
   else if (vers == calc::v4)
      echgtrn_cu2<calc::V4>();
   else if (vers == calc::v5)
      echgtrn_cu2<calc::V5>();
   else if (vers == calc::v6)
      echgtrn_cu2<calc::V6>();
}
}
