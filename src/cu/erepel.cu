#include "add.h"
#include "ff/energy.h"
#include "ff/hippo/erepel.h"
#include "ff/image.h"
#include "ff/switch.h"
#include "launch.h"
#include "mod/elecamoeba.h"
#include "mod/nblist.h"
#include "mod/repel.h"
#include "seq/bsplgen.h"
#include "seq/damp_hippo.h"
#include "seq/pair_repel.h"
#include "seq/triangle.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.3
template <class Ver>
__global__
void erepel_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nr, energy_buffer restrict er,
   virial_buffer restrict vr, grad_prec* restrict gx, grad_prec* restrict gy,
   grad_prec* restrict gz, real cut, real off, const unsigned* restrict rinfo, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale, const real* restrict x,
   const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted,
   int nakpl, const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   real* restrict trqx, real* restrict trqy, real* restrict trqz, const real (*restrict rpole)[10],
   const real* restrict sizpr, const real* restrict elepr, const real* restrict dmppr)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int nrtl;
   if CONSTEXPR (do_a) {
      nrtl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec ertl;
   if CONSTEXPR (do_e) {
      ertl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vrtlxx, vrtlyx, vrtlzx, vrtlyy, vrtlzy, vrtlzz;
   if CONSTEXPR (do_v) {
      vrtlxx = 0;
      vrtlyx = 0;
      vrtlzx = 0;
      vrtlyy = 0;
      vrtlzy = 0;
      vrtlzz = 0;
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
   real gxk;
   real gyk;
   real gzk;
   real txk;
   real tyk;
   real tzk;
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
   __shared__ real sizi[BLOCK_DIM];
   __shared__ real dmpi[BLOCK_DIM];
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
   real sizk;
   real dmpk;
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
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
      }

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii];

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
      sizi[klane] = sizpr[i];
      dmpi[klane] = dmppr[i];
      vali[klane] = elepr[i];
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
      sizk = sizpr[k];
      dmpk = dmppr[k];
      valk = elepr[k];

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];

      real e;
      PairRepelGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_repel<do_g>( //
            r2, scalea, cut, off, xr, yr, zr, sizi[klane], dmpi[klane], vali[klane], ci[klane],
            dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane],
            qiyz[klane], qizz[klane], sizk, dmpk, valk, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy,
            qkyz, qkzz, e, pgrad);

         if CONSTEXPR (do_a)
            if (e != 0)
               nrtl += 1;
         if CONSTEXPR (do_e)
            ertl += floatTo<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            gxi[klane] += pgrad.frcx;
            gyi[klane] += pgrad.frcy;
            gzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;

            txi[klane] += pgrad.ttqi[0];
            tyi[klane] += pgrad.ttqi[1];
            tzi[klane] += pgrad.ttqi[2];
            txk += pgrad.ttqk[0];
            tyk += pgrad.ttqk[1];
            tzk += pgrad.ttqk[2];
         }
         if CONSTEXPR (do_v) {
            vrtlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
            vrtlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
            vrtlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
            vrtlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
            vrtlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
            vrtlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
         }
      } // end if (include)

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], trqx, i);
         atomic_add(tyi[threadIdx.x], trqy, i);
         atomic_add(tzi[threadIdx.x], trqz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
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
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
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
      sizi[threadIdx.x] = sizpr[i];
      dmpi[threadIdx.x] = dmppr[i];
      vali[threadIdx.x] = elepr[i];
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
      sizk = sizpr[k];
      dmpk = dmppr[k];
      valk = elepr[k];

      unsigned int rinfo0 = rinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (rinfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];

         real e;
         PairRepelGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_repel<do_g>( //
               r2, scalea, cut, off, xr, yr, zr, sizi[klane], dmpi[klane], vali[klane], ci[klane],
               dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
               qiyy[klane], qiyz[klane], qizz[klane], sizk, dmpk, valk, ck, dkx, dky, dkz, qkxx,
               qkxy, qkxz, qkyy, qkyz, qkzz, e, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nrtl += 1;
            if CONSTEXPR (do_e)
               ertl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               gxi[klane] += pgrad.frcx;
               gyi[klane] += pgrad.frcy;
               gzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               txi[klane] += pgrad.ttqi[0];
               tyi[klane] += pgrad.ttqi[1];
               tzi[klane] += pgrad.ttqi[2];
               txk += pgrad.ttqk[0];
               tyk += pgrad.ttqk[1];
               tzk += pgrad.ttqk[2];
            }
            if CONSTEXPR (do_v) {
               vrtlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
               vrtlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vrtlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vrtlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
               vrtlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vrtlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
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
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
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
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
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
      sizi[threadIdx.x] = sizpr[i];
      dmpi[threadIdx.x] = dmppr[i];
      vali[threadIdx.x] = elepr[i];
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
      sizk = sizpr[k];
      dmpk = dmppr[k];
      valk = elepr[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real scalea = 1;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];

         real e;
         PairRepelGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_repel<do_g>( //
               r2, scalea, cut, off, xr, yr, zr, sizi[klane], dmpi[klane], vali[klane], ci[klane],
               dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
               qiyy[klane], qiyz[klane], qizz[klane], sizk, dmpk, valk, ck, dkx, dky, dkz, qkxx,
               qkxy, qkxz, qkyy, qkyz, qkzz, e, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nrtl += 1;
            if CONSTEXPR (do_e)
               ertl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               gxi[klane] += pgrad.frcx;
               gyi[klane] += pgrad.frcy;
               gzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               txi[klane] += pgrad.ttqi[0];
               tyi[klane] += pgrad.ttqi[1];
               tzi[klane] += pgrad.ttqi[2];
               txk += pgrad.ttqk[0];
               tyk += pgrad.ttqk[1];
               tzk += pgrad.ttqk[2];
            }
            if CONSTEXPR (do_v) {
               vrtlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
               vrtlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vrtlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vrtlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
               vrtlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vrtlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
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
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(nrtl, nr, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(ertl, er, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vrtlxx, vrtlyx, vrtlzx, vrtlyy, vrtlzy, vrtlzz, vr, ithread);
   }
}

template <class Ver>
void erepel_cu2()
{
   const auto& st = *mspatial_v2_unit;
   real cut = switchCut(Switch::REPULS);
   real off = switchOff(Switch::REPULS);

   int ngrid = gpuGridSize(BLOCK_DIM);
   erepel_cu1<Ver><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nrep, er, vir_er, derx,
      dery, derz, cut, off, st.si2.bit0, nrepexclude, repexclude, repexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, trqx, trqy, trqz, rpole, sizpr,
      elepr, dmppr);
}

void erepel_cu(int vers)
{
   if (vers == calc::v0)
      erepel_cu2<calc::V0>();
   else if (vers == calc::v1)
      erepel_cu2<calc::V1>();
   else if (vers == calc::v3)
      erepel_cu2<calc::V3>();
   else if (vers == calc::v4)
      erepel_cu2<calc::V4>();
   else if (vers == calc::v5)
      erepel_cu2<calc::V5>();
   else if (vers == calc::v6)
      erepel_cu2<calc::V6>();
}
}
