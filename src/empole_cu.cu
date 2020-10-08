#include "add.h"
#include "empole.h"
#include "empole_self.h"
#include "glob.mplar.h"
#include "glob.pme.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "seq_pair_mpole.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
template <class Ver, class ETYP>
__global__
void empole_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nem,
                energy_buffer restrict em, virial_buffer restrict vem,
                grad_prec* restrict gx, grad_prec* restrict gy,
                grad_prec* restrict gz, real off,
                const unsigned* restrict mdpuinfo, int nexclude,
                const int (*restrict exclude)[2],
                const real (*restrict exclude_scale)[4], const real* restrict x,
                const real* restrict y, const real* restrict z,
                const Spatial::SortedAtom* restrict sorted, int nakpl,
                const int* restrict iakpl, int niak, const int* restrict iak,
                const int* restrict lst, real* restrict trqx,
                real* restrict trqy, real* restrict trqz,
                const real (*restrict rpole)[10], real f, real aewald)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   int nemtl;
   if CONSTEXPR (do_a) {
      nemtl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec emtl;
   if CONSTEXPR (do_e) {
      emtl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vemtlxx, vemtlyx, vemtlzx, vemtlyy, vemtlzy, vemtlzz;
   if CONSTEXPR (do_v) {
      vemtlxx = 0;
      vemtlyx = 0;
      vemtlzx = 0;
      vemtlyy = 0;
      vemtlzy = 0;
      vemtlzz = 0;
   }


   __shared__ real shxi[BLOCK_DIM];
   __shared__ real shyi[BLOCK_DIM];
   __shared__ real shzi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real shfrcxi[BLOCK_DIM];
   __shared__ real shfrcyi[BLOCK_DIM];
   __shared__ real shfrczi[BLOCK_DIM];
   __shared__ real shtrqxi[BLOCK_DIM];
   __shared__ real shtrqyi[BLOCK_DIM];
   __shared__ real shtrqzi[BLOCK_DIM];
   real frcxk;
   real frcyk;
   real frczk;
   real trqxk;
   real trqyk;
   real trqzk;
   __shared__ real shci[BLOCK_DIM];
   __shared__ real shdix[BLOCK_DIM];
   __shared__ real shdiy[BLOCK_DIM];
   __shared__ real shdiz[BLOCK_DIM];
   __shared__ real shqixx[BLOCK_DIM];
   __shared__ real shqixy[BLOCK_DIM];
   __shared__ real shqixz[BLOCK_DIM];
   __shared__ real shqiyy[BLOCK_DIM];
   __shared__ real shqiyz[BLOCK_DIM];
   __shared__ real shqizz[BLOCK_DIM];
   real ck;
   real dkx;
   real dky;
   real dkz;
   real qkxx;
   real qkxy;
   real qkxz;
   real qkyy;
   real qkyz;
   real qkzz;


   //* /
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         shfrcxi[threadIdx.x] = 0;
         shfrcyi[threadIdx.x] = 0;
         shfrczi[threadIdx.x] = 0;
         shtrqxi[threadIdx.x] = 0;
         shtrqyi[threadIdx.x] = 0;
         shtrqzi[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         trqxk = 0;
         trqyk = 0;
         trqzk = 0;
      }


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii][0];


      real xi = x[shi];
      real yi = y[shi];
      real zi = z[shi];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      real ci = rpole[shi][mpl_pme_0];
      real dix = rpole[shi][mpl_pme_x];
      real diy = rpole[shi][mpl_pme_y];
      real diz = rpole[shi][mpl_pme_z];
      real qixx = rpole[shi][mpl_pme_xx];
      real qixy = rpole[shi][mpl_pme_xy];
      real qixz = rpole[shi][mpl_pme_xz];
      real qiyy = rpole[shi][mpl_pme_yy];
      real qiyz = rpole[shi][mpl_pme_yz];
      real qizz = rpole[shi][mpl_pme_zz];
      ck = rpole[k][mpl_pme_0];
      dkx = rpole[k][mpl_pme_x];
      dky = rpole[k][mpl_pme_y];
      dkz = rpole[k][mpl_pme_z];
      qkxx = rpole[k][mpl_pme_xx];
      qkxy = rpole[k][mpl_pme_xy];
      qkxz = rpole[k][mpl_pme_xz];
      qkyy = rpole[k][mpl_pme_yy];
      qkyz = rpole[k][mpl_pme_yz];
      qkzz = rpole[k][mpl_pme_zz];


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real e, vxx, vyx, vzx, vyy, vzy, vzz;
         real e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1;
         pair_mpole_v2<Ver, ETYP>(
            r2, xr, yr, zr, 1, ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz,
            qizz, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, f,
            aewald, shfrcxi[klane], shfrcyi[klane], shfrczi[klane], frcxk,
            frcyk, frczk, shtrqxi[klane], shtrqyi[klane], shtrqzi[klane], trqxk,
            trqyk, trqzk, e, vxx, vyx, vzx, vyy, vzy, vzz);
         pair_mpole_v2<Ver, NON_EWALD>(
            r2, xr, yr, zr, scalea - 1, ci, dix, diy, diz, qixx, qixy, qixz,
            qiyy, qiyz, qizz, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz,
            qkzz, f, aewald, shfrcxi[klane], shfrcyi[klane], shfrczi[klane],
            frcxk, frcyk, frczk, shtrqxi[klane], shtrqyi[klane], shtrqzi[klane],
            trqxk, trqyk, trqzk, e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1);
         if CONSTEXPR (do_e) {
            e = e + e1;
            emtl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalea != 0 and e != 0)
                  nemtl += 1;
            }
         }
         if CONSTEXPR (do_v) {
            vemtlxx += cvt_to<vbuf_prec>(vxx + vxx1);
            vemtlyx += cvt_to<vbuf_prec>(vyx + vyx1);
            vemtlzx += cvt_to<vbuf_prec>(vzx + vzx1);
            vemtlyy += cvt_to<vbuf_prec>(vyy + vyy1);
            vemtlzy += cvt_to<vbuf_prec>(vzy + vzy1);
            vemtlzz += cvt_to<vbuf_prec>(vzz + vzz1);
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shfrcxi[threadIdx.x], gx, shi);
         atomic_add(shfrcyi[threadIdx.x], gy, shi);
         atomic_add(shfrczi[threadIdx.x], gz, shi);
         atomic_add(shtrqxi[threadIdx.x], trqx, shi);
         atomic_add(shtrqyi[threadIdx.x], trqy, shi);
         atomic_add(shtrqzi[threadIdx.x], trqz, shi);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(trqxk, trqx, k);
         atomic_add(trqyk, trqy, k);
         atomic_add(trqzk, trqz, k);
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shfrcxi[threadIdx.x] = 0;
         shfrcyi[threadIdx.x] = 0;
         shfrczi[threadIdx.x] = 0;
         shtrqxi[threadIdx.x] = 0;
         shtrqyi[threadIdx.x] = 0;
         shtrqzi[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         trqxk = 0;
         trqyk = 0;
         trqzk = 0;
      }


      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);


      int shiid = ty * WARP_SIZE + ilane;
      int shatomi = min(shiid, n - 1);
      int shi = sorted[shatomi].unsorted;
      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      int k = sorted[atomk].unsorted;
      shxi[threadIdx.x] = sorted[shatomi].x;
      shyi[threadIdx.x] = sorted[shatomi].y;
      shzi[threadIdx.x] = sorted[shatomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;


      shci[threadIdx.x] = rpole[shi][mpl_pme_0];
      shdix[threadIdx.x] = rpole[shi][mpl_pme_x];
      shdiy[threadIdx.x] = rpole[shi][mpl_pme_y];
      shdiz[threadIdx.x] = rpole[shi][mpl_pme_z];
      shqixx[threadIdx.x] = rpole[shi][mpl_pme_xx];
      shqixy[threadIdx.x] = rpole[shi][mpl_pme_xy];
      shqixz[threadIdx.x] = rpole[shi][mpl_pme_xz];
      shqiyy[threadIdx.x] = rpole[shi][mpl_pme_yy];
      shqiyz[threadIdx.x] = rpole[shi][mpl_pme_yz];
      shqizz[threadIdx.x] = rpole[shi][mpl_pme_zz];
      ck = rpole[k][mpl_pme_0];
      dkx = rpole[k][mpl_pme_x];
      dky = rpole[k][mpl_pme_y];
      dkz = rpole[k][mpl_pme_z];
      qkxx = rpole[k][mpl_pme_xx];
      qkxy = rpole[k][mpl_pme_xy];
      qkxz = rpole[k][mpl_pme_xz];
      qkyy = rpole[k][mpl_pme_yy];
      qkyz = rpole[k][mpl_pme_yz];
      qkzz = rpole[k][mpl_pme_zz];


      unsigned int mdpuinfo0 = mdpuinfo[iw * WARP_SIZE + ilane];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int srcmask = 1 << srclane;
         int klane = srclane + threadIdx.x - ilane;
         int iid = shiid;
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real ci = shci[klane];
         real dix = shdix[klane];
         real diy = shdiy[klane];
         real diz = shdiz[klane];
         real qixx = shqixx[klane];
         real qixy = shqixy[klane];
         real qixz = shqixz[klane];
         real qiyy = shqiyy[klane];
         real qiyz = shqiyz[klane];
         real qizz = shqizz[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (mdpuinfo0 & srcmask) == 0;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
            pair_mpole_v2<Ver, ETYP>(
               r2, xr, yr, zr, 1, ci, dix, diy, diz, qixx, qixy, qixz, qiyy,
               qiyz, qizz, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz,
               qkzz, f, aewald, shfrcxi[klane], shfrcyi[klane], shfrczi[klane],
               frcxk, frcyk, frczk, shtrqxi[klane], shtrqyi[klane],
               shtrqzi[klane], trqxk, trqyk, trqzk, e, vxx, vyx, vzx, vyy, vzy,
               vzz);
            if CONSTEXPR (do_e) {
               emtl += cvt_to<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nemtl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               vemtlxx += cvt_to<vbuf_prec>(vxx);
               vemtlyx += cvt_to<vbuf_prec>(vyx);
               vemtlzx += cvt_to<vbuf_prec>(vzx);
               vemtlyy += cvt_to<vbuf_prec>(vyy);
               vemtlzy += cvt_to<vbuf_prec>(vzy);
               vemtlzz += cvt_to<vbuf_prec>(vzz);
            }
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(shfrcxi[threadIdx.x], gx, shi);
         atomic_add(shfrcyi[threadIdx.x], gy, shi);
         atomic_add(shfrczi[threadIdx.x], gz, shi);
         atomic_add(shtrqxi[threadIdx.x], trqx, shi);
         atomic_add(shtrqyi[threadIdx.x], trqy, shi);
         atomic_add(shtrqzi[threadIdx.x], trqz, shi);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(trqxk, trqx, k);
         atomic_add(trqyk, trqy, k);
         atomic_add(trqzk, trqz, k);
      }
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shfrcxi[threadIdx.x] = 0;
         shfrcyi[threadIdx.x] = 0;
         shfrczi[threadIdx.x] = 0;
         shtrqxi[threadIdx.x] = 0;
         shtrqyi[threadIdx.x] = 0;
         shtrqzi[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         trqxk = 0;
         trqyk = 0;
         trqzk = 0;
      }


      int ty = iak[iw];
      int shatomi = ty * WARP_SIZE + ilane;
      int shi = sorted[shatomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      shxi[threadIdx.x] = sorted[shatomi].x;
      shyi[threadIdx.x] = sorted[shatomi].y;
      shzi[threadIdx.x] = sorted[shatomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;


      shci[threadIdx.x] = rpole[shi][mpl_pme_0];
      shdix[threadIdx.x] = rpole[shi][mpl_pme_x];
      shdiy[threadIdx.x] = rpole[shi][mpl_pme_y];
      shdiz[threadIdx.x] = rpole[shi][mpl_pme_z];
      shqixx[threadIdx.x] = rpole[shi][mpl_pme_xx];
      shqixy[threadIdx.x] = rpole[shi][mpl_pme_xy];
      shqixz[threadIdx.x] = rpole[shi][mpl_pme_xz];
      shqiyy[threadIdx.x] = rpole[shi][mpl_pme_yy];
      shqiyz[threadIdx.x] = rpole[shi][mpl_pme_yz];
      shqizz[threadIdx.x] = rpole[shi][mpl_pme_zz];
      ck = rpole[k][mpl_pme_0];
      dkx = rpole[k][mpl_pme_x];
      dky = rpole[k][mpl_pme_y];
      dkz = rpole[k][mpl_pme_z];
      qkxx = rpole[k][mpl_pme_xx];
      qkxy = rpole[k][mpl_pme_xy];
      qkxz = rpole[k][mpl_pme_xz];
      qkyy = rpole[k][mpl_pme_yy];
      qkyz = rpole[k][mpl_pme_yz];
      qkzz = rpole[k][mpl_pme_zz];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real ci = shci[klane];
         real dix = shdix[klane];
         real diy = shdiy[klane];
         real diz = shdiz[klane];
         real qixx = shqixx[klane];
         real qixy = shqixy[klane];
         real qixz = shqixz[klane];
         real qiyy = shqiyy[klane];
         real qiyz = shqiyz[klane];
         real qizz = shqizz[klane];


         bool incl = atomk > 0;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
            pair_mpole_v2<Ver, ETYP>(
               r2, xr, yr, zr, 1, ci, dix, diy, diz, qixx, qixy, qixz, qiyy,
               qiyz, qizz, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz,
               qkzz, f, aewald, shfrcxi[klane], shfrcyi[klane], shfrczi[klane],
               frcxk, frcyk, frczk, shtrqxi[klane], shtrqyi[klane],
               shtrqzi[klane], trqxk, trqyk, trqzk, e, vxx, vyx, vzx, vyy, vzy,
               vzz);
            if CONSTEXPR (do_e) {
               emtl += cvt_to<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nemtl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               vemtlxx += cvt_to<vbuf_prec>(vxx);
               vemtlyx += cvt_to<vbuf_prec>(vyx);
               vemtlzx += cvt_to<vbuf_prec>(vzx);
               vemtlyy += cvt_to<vbuf_prec>(vyy);
               vemtlzy += cvt_to<vbuf_prec>(vzy);
               vemtlzz += cvt_to<vbuf_prec>(vzz);
            }
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(shfrcxi[threadIdx.x], gx, shi);
         atomic_add(shfrcyi[threadIdx.x], gy, shi);
         atomic_add(shfrczi[threadIdx.x], gz, shi);
         atomic_add(shtrqxi[threadIdx.x], trqx, shi);
         atomic_add(shtrqyi[threadIdx.x], trqy, shi);
         atomic_add(shtrqzi[threadIdx.x], trqz, shi);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(trqxk, trqx, k);
         atomic_add(trqyk, trqy, k);
         atomic_add(trqzk, trqz, k);
      }
   }
   // */


   if CONSTEXPR (do_a) {
      atomic_add(nemtl, nem, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(emtl, em, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vemtlxx, vemtlyx, vemtlzx, vemtlyy, vemtlzy, vemtlzz, vem,
                 ithread);
   }
} // generated by ComplexKernelBuilder (ck.py) 1.5.2


template <class Ver, class ETYP>
void empole_cu()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;


   const auto& st = *mspatial_v2_unit;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;


      if CONSTEXPR (do_e) {
         launch_k1s(nonblk, n, empole_self_cu<do_a>, //
                    nem, em, rpole, n, f, aewald);
      }
   }
   int ngrid = get_grid_size(BLOCK_DIM);
   empole_cu1<Ver, ETYP><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, nem, em, vir_em, demx, demy, demz, off,
      st.si1.bit0, nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, trqx, trqy,
      trqz, rpole, f, aewald);
}


void empole_nonewald_cu(int vers)
{
   if (vers == calc::v0) {
      empole_cu<calc::V0, NON_EWALD>();
   } else if (vers == calc::v1) {
      empole_cu<calc::V1, NON_EWALD>();
   } else if (vers == calc::v3) {
      empole_cu<calc::V3, NON_EWALD>();
   } else if (vers == calc::v4) {
      empole_cu<calc::V4, NON_EWALD>();
   } else if (vers == calc::v5) {
      empole_cu<calc::V5, NON_EWALD>();
   } else if (vers == calc::v6) {
      empole_cu<calc::V6, NON_EWALD>();
   }
}


void empole_ewald_real_self_cu(int vers)
{
   if (vers == calc::v0) {
      empole_cu<calc::V0, EWALD>();
   } else if (vers == calc::v1) {
      empole_cu<calc::V1, EWALD>();
   } else if (vers == calc::v3) {
      empole_cu<calc::V3, EWALD>();
   } else if (vers == calc::v4) {
      empole_cu<calc::V4, EWALD>();
   } else if (vers == calc::v5) {
      empole_cu<calc::V5, EWALD>();
   } else if (vers == calc::v6) {
      empole_cu<calc::V6, EWALD>();
   }
}
}
