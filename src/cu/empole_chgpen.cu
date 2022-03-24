#include "add.h"
#include "empole_chgpen_self.h"
#include "ff/hippo/cflux.h"
#include "ff/hippo/empole_chgpen.h"
#include "ff/image.h"
#include "ff/pmestuf.h"
#include "ff/switch.h"
#include "glob/pme.h"
#include "glob/spatial.h"
#include "launch.h"
#include "md/inc.h"
#include "mod/elecpchg.h"
#include "seq/bsplgen.h"
#include "seq/pair_mpole_chgpen.h"
#include "seq/triangle.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.3
template <class Ver, class ETYP, bool CFLX>
__global__
void empole_chgpen_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nem,
   energy_buffer restrict em, virial_buffer restrict vem, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, real off, const unsigned* restrict minfo,
   int nexclude, const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, real* restrict trqx, real* restrict trqy,
   real* restrict trqz, real* restrict pot, const real (*restrict rpole)[10], real* restrict pcore,
   real* restrict pval, const real* restrict palpha, real aewald, real f)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;
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
   __shared__ real xi[BLOCK_DIM];
   __shared__ real yi[BLOCK_DIM];
   __shared__ real zi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real gxi[BLOCK_DIM];
   __shared__ real gyi[BLOCK_DIM];
   __shared__ real gzi[BLOCK_DIM];
   __shared__ real txi[BLOCK_DIM];
   __shared__ real tyi[BLOCK_DIM];
   __shared__ real tzi[BLOCK_DIM];
   __shared__ real poti[BLOCK_DIM];
   real gxk;
   real gyk;
   real gzk;
   real txk;
   real tyk;
   real tzk;
   real potk;
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
   __shared__ real corei[BLOCK_DIM];
   __shared__ real alphai[BLOCK_DIM];
   __shared__ real vali[BLOCK_DIM];
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
   real corek;
   real alphak;
   real valk;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      if CONSTEXPR (do_g) {
         gxi[threadIdx.x] = 0;
         gyi[threadIdx.x] = 0;
         gzi[threadIdx.x] = 0;
         txi[threadIdx.x] = 0;
         tyi[threadIdx.x] = 0;
         tzi[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
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
      corei[klane] = pcore[i];
      alphai[klane] = palpha[i];
      vali[klane] = pval[i];
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];

      real e;
      real pota, potb;
      PairMPoleGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_mpole_chgpen<do_e, do_g, ETYP, CFLX>(r2, xr, yr, zr, scalea, //
            ci[klane], dix[klane], diy[klane], diz[klane], corei[klane], vali[klane],
            alphai[klane], //
            qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
            qizz[klane],                            //
            ck, dkx, dky, dkz, corek, valk, alphak, //
            qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,     //
            f, aewald, e, pota, potb, pgrad);

         if CONSTEXPR (do_a)
            if (e != 0 and scalea != 0)
               nemtl += 1;
         if CONSTEXPR (do_e)
            emtl += floatTo<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            gxi[klane] += pgrad.frcx;
            gyi[klane] += pgrad.frcy;
            gzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;

            txi[klane] += pgrad.ttmi[0];
            tyi[klane] += pgrad.ttmi[1];
            tzi[klane] += pgrad.ttmi[2];
            txk += pgrad.ttmk[0];
            tyk += pgrad.ttmk[1];
            tzk += pgrad.ttmk[2];
         }
         if CONSTEXPR (do_v) {
            vemtlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
            vemtlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
            vemtlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
            vemtlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
            vemtlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
            vemtlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
         }
         if CONSTEXPR (CFLX) {
            poti[klane] += pota;
            potk += potb;
         }
      } // end if (include)

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], trqx, i);
         atomic_add(tyi[threadIdx.x], trqy, i);
         atomic_add(tzi[threadIdx.x], trqz, i);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi[threadIdx.x] = 0;
         gyi[threadIdx.x] = 0;
         gzi[threadIdx.x] = 0;
         txi[threadIdx.x] = 0;
         tyi[threadIdx.x] = 0;
         tzi[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
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
      corei[threadIdx.x] = pcore[i];
      alphai[threadIdx.x] = palpha[i];
      vali[threadIdx.x] = pval[i];
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];

      unsigned int minfo0 = minfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (minfo0 & srcmask) == 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];

         real e;
         real pota, potb;
         PairMPoleGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_mpole_chgpen<do_e, do_g, ETYP, CFLX>(r2, xr, yr, zr, 1, //
               ci[klane], dix[klane], diy[klane], diz[klane], corei[klane], vali[klane],
               alphai[klane], //
               qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
               qizz[klane],                            //
               ck, dkx, dky, dkz, corek, valk, alphak, //
               qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,     //
               f, aewald, e, pota, potb, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nemtl += 1;
            if CONSTEXPR (do_e)
               emtl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               gxi[klane] += pgrad.frcx;
               gyi[klane] += pgrad.frcy;
               gzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               txi[klane] += pgrad.ttmi[0];
               tyi[klane] += pgrad.ttmi[1];
               tzi[klane] += pgrad.ttmi[2];
               txk += pgrad.ttmk[0];
               tyk += pgrad.ttmk[1];
               tzk += pgrad.ttmk[2];
            }
            if CONSTEXPR (do_v) {
               vemtlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
               vemtlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vemtlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vemtlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
               vemtlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vemtlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               poti[klane] += pota;
               potk += potb;
            }
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], trqx, i);
         atomic_add(tyi[threadIdx.x], trqy, i);
         atomic_add(tzi[threadIdx.x], trqz, i);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi[threadIdx.x] = 0;
         gyi[threadIdx.x] = 0;
         gzi[threadIdx.x] = 0;
         txi[threadIdx.x] = 0;
         tyi[threadIdx.x] = 0;
         tzi[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
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
      corei[threadIdx.x] = pcore[i];
      alphai[threadIdx.x] = palpha[i];
      vali[threadIdx.x] = pval[i];
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];

         real e;
         real pota, potb;
         PairMPoleGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_mpole_chgpen<do_e, do_g, ETYP, CFLX>(r2, xr, yr, zr, 1, //
               ci[klane], dix[klane], diy[klane], diz[klane], corei[klane], vali[klane],
               alphai[klane], //
               qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
               qizz[klane],                            //
               ck, dkx, dky, dkz, corek, valk, alphak, //
               qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,     //
               f, aewald, e, pota, potb, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nemtl += 1;
            if CONSTEXPR (do_e)
               emtl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               gxi[klane] += pgrad.frcx;
               gyi[klane] += pgrad.frcy;
               gzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               txi[klane] += pgrad.ttmi[0];
               tyi[klane] += pgrad.ttmi[1];
               tzi[klane] += pgrad.ttmi[2];
               txk += pgrad.ttmk[0];
               tyk += pgrad.ttmk[1];
               tzk += pgrad.ttmk[2];
            }
            if CONSTEXPR (do_v) {
               vemtlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
               vemtlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vemtlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vemtlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
               vemtlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vemtlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               poti[klane] += pota;
               potk += potb;
            }
         } // end if (include)
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], trqx, i);
         atomic_add(tyi[threadIdx.x], trqy, i);
         atomic_add(tzi[threadIdx.x], trqz, i);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
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

template <class Ver, class ETYP, int CFLX>
void empole_chgpen_cu()
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
      launch_k1b(g::s0, n, empole_chgpen_self_cu<do_a, do_e, CFLX>, //
         nem, em, rpole, pot, n, f, aewald);
   }

   int ngrid = gpuGridSize(BLOCK_DIM);
   empole_chgpen_cu1<Ver, ETYP, CFLX><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nem,
      em, vir_em, demx, demy, demz, off, st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale,
      st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, trqx, trqy, trqz,
      pot, rpole, pcore, pval, palpha, aewald, f);
}

void empole_chgpen_nonewald_cu(int vers, int use_cf)
{
   if (use_cf) {
      if (vers == calc::v0) {
         // empole_chgpen_cu<calc::V0, NON_EWALD, 1>();
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, NON_EWALD, 1>();
      } else if (vers == calc::v3) {
         // empole_chgpen_cu<calc::V3, NON_EWALD, 1>();
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, NON_EWALD, 1>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, NON_EWALD, 1>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, NON_EWALD, 1>();
      }
   } else {
      if (vers == calc::v0) {
         empole_chgpen_cu<calc::V0, NON_EWALD, 0>();
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, NON_EWALD, 0>();
      } else if (vers == calc::v3) {
         empole_chgpen_cu<calc::V3, NON_EWALD, 0>();
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, NON_EWALD, 0>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, NON_EWALD, 0>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, NON_EWALD, 0>();
      }
   }
}

void empole_chgpen_ewald_real_self_cu(int vers, int use_cf)
{
   if (use_cf) {
      if (vers == calc::v0) {
         // empole_chgpen_cu<calc::V0, EWALD, 1>();
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, EWALD, 1>();
      } else if (vers == calc::v3) {
         // empole_chgpen_cu<calc::V3, EWALD, 1>();
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, EWALD, 1>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, EWALD, 1>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, EWALD, 1>();
      }
   } else {
      if (vers == calc::v0) {
         empole_chgpen_cu<calc::V0, EWALD, 0>();
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, EWALD, 0>();
      } else if (vers == calc::v3) {
         empole_chgpen_cu<calc::V3, EWALD, 0>();
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, EWALD, 0>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, EWALD, 0>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, EWALD, 0>();
      }
   }
}
}
