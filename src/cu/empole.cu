#include "add.h"
#include "empole_self.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/amoeba/empole.h"
#include "ff/energy.h"
#include "ff/image.h"
#include "ff/pchg/echarge.h"
#include "ff/pme.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "launch.h"
#include "seq/pair_mpole.h"
#include "seq/triangle.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.2
template <class Ver, class ETYP>
__global__
void empole_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nem, EnergyBuffer restrict em,
   VirialBuffer restrict vem, grad_prec* restrict gx, grad_prec* restrict gy,
   grad_prec* restrict gz, real off, const unsigned* restrict mdpuinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[4],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, real* restrict trqx, real* restrict trqy,
   real* restrict trqz, const real (*restrict rpole)[10], real f, real aewald)
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
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec emtl;
   if CONSTEXPR (do_e) {
      emtl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec vemtlxx, vemtlyx, vemtlzx, vemtlyy, vemtlzy, vemtlzz;
   if CONSTEXPR (do_v) {
      vemtlxx = 0;
      vemtlyx = 0;
      vemtlzx = 0;
      vemtlyy = 0;
      vemtlzy = 0;
      vemtlzz = 0;
   }
   __shared__ real xi[BLOCK_DIM];
   __shared__ real yi[BLOCK_DIM];
   __shared__ real zi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real frcxi[BLOCK_DIM];
   __shared__ real frcyi[BLOCK_DIM];
   __shared__ real frczi[BLOCK_DIM];
   __shared__ real trqxi[BLOCK_DIM];
   __shared__ real trqyi[BLOCK_DIM];
   __shared__ real trqzi[BLOCK_DIM];
   real frcxk;
   real frcyk;
   real frczk;
   real trqxk;
   real trqyk;
   real trqzk;
   __shared__ real ci[BLOCK_DIM];
   __shared__ real dix[BLOCK_DIM];
   __shared__ real diy[BLOCK_DIM];
   __shared__ real diz[BLOCK_DIM];
   __shared__ real qixx[BLOCK_DIM];
   __shared__ real qixy[BLOCK_DIM];
   __shared__ real qixz[BLOCK_DIM];
   __shared__ real qiyy[BLOCK_DIM];
   __shared__ real qiyz[BLOCK_DIM];
   __shared__ real qizz[BLOCK_DIM];
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
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      if CONSTEXPR (do_g) {
         frcxi[threadIdx.x] = 0;
         frcyi[threadIdx.x] = 0;
         frczi[threadIdx.x] = 0;
         trqxi[threadIdx.x] = 0;
         trqyi[threadIdx.x] = 0;
         trqzi[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         trqxk = 0;
         trqyk = 0;
         trqzk = 0;
      }

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii][0];

      xi[klane] = x[i];
      yi[klane] = y[i];
      zi[klane] = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      ci[klane] = rpole[i][mpl_pme_0];
      dix[klane] = rpole[i][mpl_pme_x];
      diy[klane] = rpole[i][mpl_pme_y];
      diz[klane] = rpole[i][mpl_pme_z];
      qixx[klane] = rpole[i][mpl_pme_xx];
      qixy[klane] = rpole[i][mpl_pme_xy];
      qixz[klane] = rpole[i][mpl_pme_xz];
      qiyy[klane] = rpole[i][mpl_pme_yy];
      qiyz[klane] = rpole[i][mpl_pme_yz];
      qizz[klane] = rpole[i][mpl_pme_zz];
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

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real e, vxx, vyx, vzx, vyy, vzy, vzz;
         real e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1;
         pair_mpole_v2<Ver, ETYP>(r2, xr, yr, zr, 1, ci[klane], dix[klane], diy[klane], diz[klane],
            qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], ck, dkx,
            dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, f, aewald, frcxi[klane], frcyi[klane],
            frczi[klane], frcxk, frcyk, frczk, trqxi[klane], trqyi[klane], trqzi[klane], trqxk,
            trqyk, trqzk, e, vxx, vyx, vzx, vyy, vzy, vzz);
         pair_mpole_v2<Ver, NON_EWALD>(r2, xr, yr, zr, scalea - 1, ci[klane], dix[klane],
            diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
            qizz[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, f, aewald,
            frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, trqxi[klane],
            trqyi[klane], trqzi[klane], trqxk, trqyk, trqzk, e1, vxx1, vyx1, vzx1, vyy1, vzy1,
            vzz1);
         if CONSTEXPR (do_e) {
            e = e + e1;
            emtl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalea != 0 and e != 0)
                  nemtl += 1;
            }
         }
         if CONSTEXPR (do_v) {
            vemtlxx += floatTo<vbuf_prec>(vxx + vxx1);
            vemtlyx += floatTo<vbuf_prec>(vyx + vyx1);
            vemtlzx += floatTo<vbuf_prec>(vzx + vzx1);
            vemtlyy += floatTo<vbuf_prec>(vyy + vyy1);
            vemtlzy += floatTo<vbuf_prec>(vzy + vzy1);
            vemtlzz += floatTo<vbuf_prec>(vzz + vzz1);
         }
      } // end if (include)

      if CONSTEXPR (do_g) {
         atomic_add(frcxi[threadIdx.x], gx, i);
         atomic_add(frcyi[threadIdx.x], gy, i);
         atomic_add(frczi[threadIdx.x], gz, i);
         atomic_add(trqxi[threadIdx.x], trqx, i);
         atomic_add(trqyi[threadIdx.x], trqy, i);
         atomic_add(trqzi[threadIdx.x], trqz, i);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(trqxk, trqx, k);
         atomic_add(trqyk, trqy, k);
         atomic_add(trqzk, trqz, k);
      }
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         frcxi[threadIdx.x] = 0;
         frcyi[threadIdx.x] = 0;
         frczi[threadIdx.x] = 0;
         trqxi[threadIdx.x] = 0;
         trqyi[threadIdx.x] = 0;
         trqzi[threadIdx.x] = 0;
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

      int iid = ty * WARP_SIZE + ilane;
      int atomi = min(iid, n - 1);
      int i = sorted[atomi].unsorted;
      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      int k = sorted[atomk].unsorted;
      xi[threadIdx.x] = sorted[atomi].x;
      yi[threadIdx.x] = sorted[atomi].y;
      zi[threadIdx.x] = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      ci[threadIdx.x] = rpole[i][mpl_pme_0];
      dix[threadIdx.x] = rpole[i][mpl_pme_x];
      diy[threadIdx.x] = rpole[i][mpl_pme_y];
      diz[threadIdx.x] = rpole[i][mpl_pme_z];
      qixx[threadIdx.x] = rpole[i][mpl_pme_xx];
      qixy[threadIdx.x] = rpole[i][mpl_pme_xy];
      qixz[threadIdx.x] = rpole[i][mpl_pme_xz];
      qiyy[threadIdx.x] = rpole[i][mpl_pme_yy];
      qiyz[threadIdx.x] = rpole[i][mpl_pme_yz];
      qizz[threadIdx.x] = rpole[i][mpl_pme_zz];
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
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (mdpuinfo0 & srcmask) == 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
            pair_mpole_v2<Ver, ETYP>(r2, xr, yr, zr, 1, ci[klane], dix[klane], diy[klane],
               diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
               qizz[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, f, aewald,
               frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, trqxi[klane],
               trqyi[klane], trqzi[klane], trqxk, trqyk, trqzk, e, vxx, vyx, vzx, vyy, vzy, vzz);
            if CONSTEXPR (do_e) {
               emtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nemtl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               vemtlxx += floatTo<vbuf_prec>(vxx);
               vemtlyx += floatTo<vbuf_prec>(vyx);
               vemtlzx += floatTo<vbuf_prec>(vzx);
               vemtlyy += floatTo<vbuf_prec>(vyy);
               vemtlzy += floatTo<vbuf_prec>(vzy);
               vemtlzz += floatTo<vbuf_prec>(vzz);
            }
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
      }

      if CONSTEXPR (do_g) {
         atomic_add(frcxi[threadIdx.x], gx, i);
         atomic_add(frcyi[threadIdx.x], gy, i);
         atomic_add(frczi[threadIdx.x], gz, i);
         atomic_add(trqxi[threadIdx.x], trqx, i);
         atomic_add(trqyi[threadIdx.x], trqy, i);
         atomic_add(trqzi[threadIdx.x], trqz, i);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(trqxk, trqx, k);
         atomic_add(trqyk, trqy, k);
         atomic_add(trqzk, trqz, k);
      }
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         frcxi[threadIdx.x] = 0;
         frcyi[threadIdx.x] = 0;
         frczi[threadIdx.x] = 0;
         trqxi[threadIdx.x] = 0;
         trqyi[threadIdx.x] = 0;
         trqzi[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         trqxk = 0;
         trqyk = 0;
         trqzk = 0;
      }

      int ty = iak[iw];
      int atomi = ty * WARP_SIZE + ilane;
      int i = sorted[atomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      xi[threadIdx.x] = sorted[atomi].x;
      yi[threadIdx.x] = sorted[atomi].y;
      zi[threadIdx.x] = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      ci[threadIdx.x] = rpole[i][mpl_pme_0];
      dix[threadIdx.x] = rpole[i][mpl_pme_x];
      diy[threadIdx.x] = rpole[i][mpl_pme_y];
      diz[threadIdx.x] = rpole[i][mpl_pme_z];
      qixx[threadIdx.x] = rpole[i][mpl_pme_xx];
      qixy[threadIdx.x] = rpole[i][mpl_pme_xy];
      qixz[threadIdx.x] = rpole[i][mpl_pme_xz];
      qiyy[threadIdx.x] = rpole[i][mpl_pme_yy];
      qiyz[threadIdx.x] = rpole[i][mpl_pme_yz];
      qizz[threadIdx.x] = rpole[i][mpl_pme_zz];
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
         bool incl = atomk > 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
            pair_mpole_v2<Ver, ETYP>(r2, xr, yr, zr, 1, ci[klane], dix[klane], diy[klane],
               diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
               qizz[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, f, aewald,
               frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, trqxi[klane],
               trqyi[klane], trqzi[klane], trqxk, trqyk, trqzk, e, vxx, vyx, vzx, vyy, vzy, vzz);
            if CONSTEXPR (do_e) {
               emtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nemtl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               vemtlxx += floatTo<vbuf_prec>(vxx);
               vemtlyx += floatTo<vbuf_prec>(vyx);
               vemtlzx += floatTo<vbuf_prec>(vzx);
               vemtlyy += floatTo<vbuf_prec>(vyy);
               vemtlzy += floatTo<vbuf_prec>(vzy);
               vemtlzz += floatTo<vbuf_prec>(vzz);
            }
         } // end if (include)
      }

      if CONSTEXPR (do_g) {
         atomic_add(frcxi[threadIdx.x], gx, i);
         atomic_add(frcyi[threadIdx.x], gy, i);
         atomic_add(frczi[threadIdx.x], gz, i);
         atomic_add(trqxi[threadIdx.x], trqx, i);
         atomic_add(trqyi[threadIdx.x], trqy, i);
         atomic_add(trqzi[threadIdx.x], trqz, i);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(trqxk, trqx, k);
         atomic_add(trqyk, trqy, k);
         atomic_add(trqzk, trqz, k);
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(nemtl, nem, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(emtl, em, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vemtlxx, vemtlyx, vemtlzx, vemtlyy, vemtlzy, vemtlzz, vem, ithread);
   }
}

template <class Ver, class ETYP>
void empole_cu()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;

   const auto& st = *mspatial_v2_unit;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switchOff(Switch::EWALD);
   else
      off = switchOff(Switch::MPOLE);

   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;

      if CONSTEXPR (do_e) {
         launch_k1b(g::s0, n, empole_self_cu<do_a>, //
            nem, em, rpole, n, f, aewald);
      }
   }
   int ngrid = gpuGridSize(BLOCK_DIM);
   empole_cu1<Ver, ETYP><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nem, em, vir_em,
      demx, demy, demz, off, st.si1.bit0, nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, trqx, trqy, trqz, rpole, f,
      aewald);
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
