#include "add.h"
#include "erepel.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "seq_bsplgen.h"
#include "seq_damprep.h"
#include "seq_pair_repel.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"

namespace tinker {
template <class Ver>
__global__
void erepel_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nr,
                energy_buffer restrict er, virial_buffer restrict vr,
                grad_prec* restrict gx, grad_prec* restrict gy,
                grad_prec* restrict gz, real cut, real off,
                const unsigned* restrict rinfo, int nexclude,
                const int (*restrict exclude)[2],
                const real* restrict exclude_scale, const real* restrict x,
                const real* restrict y, const real* restrict z,
                const Spatial::SortedAtom* restrict sorted, int nakpl,
                const int* restrict iakpl, int niak, const int* restrict iak,
                const int* restrict lst, real* restrict trqx,
                real* restrict trqy, real* restrict trqz,
                const real (*restrict rpole)[10], const real* restrict sizpr,
                const real* restrict elepr, const real* restrict dmppr)
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


   __shared__ real shxi[BLOCK_DIM];
   __shared__ real shyi[BLOCK_DIM];
   __shared__ real shzi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real shgxi[BLOCK_DIM];
   __shared__ real shgyi[BLOCK_DIM];
   __shared__ real shgzi[BLOCK_DIM];
   __shared__ real shtxi[BLOCK_DIM];
   __shared__ real shtyi[BLOCK_DIM];
   __shared__ real shtzi[BLOCK_DIM];
   real gxk;
   real gyk;
   real gzk;
   real txk;
   real tyk;
   real tzk;
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
   __shared__ real shsizi[BLOCK_DIM];
   __shared__ real shdmpi[BLOCK_DIM];
   __shared__ real shvali[BLOCK_DIM];
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
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
      }


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii];


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
      real sizi = sizpr[shi];
      real dmpi = dmppr[shi];
      real vali = elepr[shi];
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


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;

      real e;
      PairRepelGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_repel<do_g>( //
            r2, scalea, cut, off, xr, yr, zr, sizi, dmpi, vali, ci, dix, diy,
            diz, qixx, qixy, qixz, qiyy, qiyz, qizz, sizk, dmpk, valk, ck, dkx,
            dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, e, pgrad);

         if CONSTEXPR (do_a)
            if (e != 0)
               nrtl += 1;
         if CONSTEXPR (do_e)
            ertl += cvt_to<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            shgxi[klane] += pgrad.frcx;
            shgyi[klane] += pgrad.frcy;
            shgzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;


            shtxi[klane] += pgrad.ttqi[0];
            shtyi[klane] += pgrad.ttqi[1];
            shtzi[klane] += pgrad.ttqi[2];
            txk += pgrad.ttqk[0];
            tyk += pgrad.ttqk[1];
            tzk += pgrad.ttqk[2];
         }
         if CONSTEXPR (do_v) {
            vrtlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
            vrtlyx +=
               cvt_to<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
            vrtlzx +=
               cvt_to<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
            vrtlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
            vrtlzy +=
               cvt_to<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
            vrtlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], trqx, shi);
         atomic_add(shtyi[threadIdx.x], trqy, shi);
         atomic_add(shtzi[threadIdx.x], trqz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
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
      shsizi[threadIdx.x] = sizpr[shi];
      shdmpi[threadIdx.x] = dmppr[shi];
      shvali[threadIdx.x] = elepr[shi];
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
         real sizi = shsizi[klane];
         real dmpi = shdmpi[klane];
         real vali = shvali[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (rinfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         PairRepelGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_repel<do_g>( //
               r2, scalea, cut, off, xr, yr, zr, sizi, dmpi, vali, ci, dix, diy,
               diz, qixx, qixy, qixz, qiyy, qiyz, qizz, sizk, dmpk, valk, ck,
               dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, e, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nrtl += 1;
            if CONSTEXPR (do_e)
               ertl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;


               shtxi[klane] += pgrad.ttqi[0];
               shtyi[klane] += pgrad.ttqi[1];
               shtzi[klane] += pgrad.ttqi[2];
               txk += pgrad.ttqk[0];
               tyk += pgrad.ttqk[1];
               tzk += pgrad.ttqk[2];
            }
            if CONSTEXPR (do_v) {
               vrtlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
               vrtlyx += cvt_to<vbuf_prec>(-0.5f *
                                           (yr * pgrad.frcx + xr * pgrad.frcy));
               vrtlzx += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcx + xr * pgrad.frcz));
               vrtlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
               vrtlzy += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcy + yr * pgrad.frcz));
               vrtlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
            }
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], trqx, shi);
         atomic_add(shtyi[threadIdx.x], trqy, shi);
         atomic_add(shtzi[threadIdx.x], trqz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
      }
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
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
      shsizi[threadIdx.x] = sizpr[shi];
      shdmpi[threadIdx.x] = dmppr[shi];
      shvali[threadIdx.x] = elepr[shi];
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
         real sizi = shsizi[klane];
         real dmpi = shdmpi[klane];
         real vali = shvali[klane];


         bool incl = atomk > 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         PairRepelGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_repel<do_g>( //
               r2, scalea, cut, off, xr, yr, zr, sizi, dmpi, vali, ci, dix, diy,
               diz, qixx, qixy, qixz, qiyy, qiyz, qizz, sizk, dmpk, valk, ck,
               dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, e, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nrtl += 1;
            if CONSTEXPR (do_e)
               ertl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;


               shtxi[klane] += pgrad.ttqi[0];
               shtyi[klane] += pgrad.ttqi[1];
               shtzi[klane] += pgrad.ttqi[2];
               txk += pgrad.ttqk[0];
               tyk += pgrad.ttqk[1];
               tzk += pgrad.ttqk[2];
            }
            if CONSTEXPR (do_v) {
               vrtlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
               vrtlyx += cvt_to<vbuf_prec>(-0.5f *
                                           (yr * pgrad.frcx + xr * pgrad.frcy));
               vrtlzx += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcx + xr * pgrad.frcz));
               vrtlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
               vrtlzy += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcy + yr * pgrad.frcz));
               vrtlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
            }
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], trqx, shi);
         atomic_add(shtyi[threadIdx.x], trqy, shi);
         atomic_add(shtzi[threadIdx.x], trqz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
      }
   }
   // */


   if CONSTEXPR (do_a) {
      atomic_add(nrtl, nr, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(ertl, er, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vrtlxx, vrtlyx, vrtlzx, vrtlyy, vrtlzy, vrtlzz, vr, ithread);
   }
} // generated by ComplexKernelBuilder (ck.py) 1.5.3


template <class Ver>
void erepel_cu2()
{
   const auto& st = *mspatial_v2_unit;
   real cut = switch_cut(switch_repuls);
   real off = switch_off(switch_repuls);


   int ngrid = get_grid_size(BLOCK_DIM);
   erepel_cu1<Ver><<<ngrid, BLOCK_DIM, 0, g::s0>>>(
      st.n, TINKER_IMAGE_ARGS, nrep, er, vir_er, derx, dery, derz, cut, off,
      st.si2.bit0, nrepexclude, repexclude, repexclude_scale, st.x, st.y, st.z,
      st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, trqx, trqy, trqz,
      rpole, sizpr, elepr, dmppr);
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
