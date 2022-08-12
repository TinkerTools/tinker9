#include "ff/amoebamod.h"
#include "ff/elec.h"
#include "ff/hippo/expol.h"
#include "ff/hippomod.h"
#include "ff/image.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "seq/add.h"
#include "seq/launch.h"
#include "seq/pair_alterpol.h"
#include "seq/triangle.h"

namespace tinker {
// ck.py Version 2.0.3

__global__
static void alterpol_cu1(int n, TINKER_IMAGE_PARAMS, real cut, real off,
   const unsigned* restrict dinfo, int nexclude, const int (*restrict exclude)[2],
   const real (*restrict exclude_scale)[3], const real* restrict x, const real* restrict y,
   const real* restrict z, const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   real (*restrict polscale)[9], const real* restrict kpep, const real* restrict prepep,
   const real* restrict dmppep, const int* restrict lpep, ExpolScr scrtyp)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   __shared__ real xi[BLOCK_DIM];
   __shared__ real yi[BLOCK_DIM];
   __shared__ real zi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real psci00[BLOCK_DIM];
   __shared__ real psci01[BLOCK_DIM];
   __shared__ real psci02[BLOCK_DIM];
   __shared__ real psci10[BLOCK_DIM];
   __shared__ real psci11[BLOCK_DIM];
   __shared__ real psci12[BLOCK_DIM];
   __shared__ real psci20[BLOCK_DIM];
   __shared__ real psci21[BLOCK_DIM];
   __shared__ real psci22[BLOCK_DIM];
   real psck00;
   real psck01;
   real psck02;
   real psck10;
   real psck11;
   real psck12;
   real psck20;
   real psck21;
   real psck22;
   __shared__ real springi[BLOCK_DIM];
   __shared__ real sizi[BLOCK_DIM];
   __shared__ real alphai[BLOCK_DIM];
   __shared__ int epli[BLOCK_DIM];
   real springk;
   real sizk;
   real alphak;
   int eplk;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      psci00[threadIdx.x] = 0;
      psci01[threadIdx.x] = 0;
      psci02[threadIdx.x] = 0;
      psci10[threadIdx.x] = 0;
      psci11[threadIdx.x] = 0;
      psci12[threadIdx.x] = 0;
      psci20[threadIdx.x] = 0;
      psci21[threadIdx.x] = 0;
      psci22[threadIdx.x] = 0;
      psck00 = 0;
      psck01 = 0;
      psck02 = 0;
      psck10 = 0;
      psck11 = 0;
      psck12 = 0;
      psck20 = 0;
      psck21 = 0;
      psck22 = 0;

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];

      xi[klane] = x[i];
      yi[klane] = y[i];
      zi[klane] = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      springi[klane] = kpep[i];
      sizi[klane] = prepep[i];
      alphai[klane] = dmppep[i];
      epli[klane] = lpep[i];
      springk = kpep[k];
      sizk = prepep[k];
      alphak = dmppep[k];
      eplk = lpep[k];

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];
      real r2 = image2(xr, yr, zr);
      if ((eplk or epli[klane]) and r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real ks2i[3][3], ks2k[3][3];
         pair_alterpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, springi[klane], sizi[klane],
            alphai[klane], springk, sizk, alphak, ks2i, ks2k);
         psci00[klane] += ks2i[0][0];
         psci01[klane] += ks2i[0][1];
         psci02[klane] += ks2i[0][2];
         psci10[klane] += ks2i[1][0];
         psci11[klane] += ks2i[1][1];
         psci12[klane] += ks2i[1][2];
         psci20[klane] += ks2i[2][0];
         psci21[klane] += ks2i[2][1];
         psci22[klane] += ks2i[2][2];
         psck00 += ks2k[0][0];
         psck01 += ks2k[0][1];
         psck02 += ks2k[0][2];
         psck10 += ks2k[1][0];
         psck11 += ks2k[1][1];
         psck12 += ks2k[1][2];
         psck20 += ks2k[2][0];
         psck21 += ks2k[2][1];
         psck22 += ks2k[2][2];
      }

      atomic_add(psci00[threadIdx.x], &polscale[i][0]);
      atomic_add(psci01[threadIdx.x], &polscale[i][1]);
      atomic_add(psci02[threadIdx.x], &polscale[i][2]);
      atomic_add(psci10[threadIdx.x], &polscale[i][3]);
      atomic_add(psci11[threadIdx.x], &polscale[i][4]);
      atomic_add(psci12[threadIdx.x], &polscale[i][5]);
      atomic_add(psci20[threadIdx.x], &polscale[i][6]);
      atomic_add(psci21[threadIdx.x], &polscale[i][7]);
      atomic_add(psci22[threadIdx.x], &polscale[i][8]);
      atomic_add(psck00, &polscale[k][0]);
      atomic_add(psck01, &polscale[k][1]);
      atomic_add(psck02, &polscale[k][2]);
      atomic_add(psck10, &polscale[k][3]);
      atomic_add(psck11, &polscale[k][4]);
      atomic_add(psck12, &polscale[k][5]);
      atomic_add(psck20, &polscale[k][6]);
      atomic_add(psck21, &polscale[k][7]);
      atomic_add(psck22, &polscale[k][8]);
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      psci00[threadIdx.x] = 0;
      psci01[threadIdx.x] = 0;
      psci02[threadIdx.x] = 0;
      psci10[threadIdx.x] = 0;
      psci11[threadIdx.x] = 0;
      psci12[threadIdx.x] = 0;
      psci20[threadIdx.x] = 0;
      psci21[threadIdx.x] = 0;
      psci22[threadIdx.x] = 0;
      psck00 = 0;
      psck01 = 0;
      psck02 = 0;
      psck10 = 0;
      psck11 = 0;
      psck12 = 0;
      psck20 = 0;
      psck21 = 0;
      psck22 = 0;

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

      springi[threadIdx.x] = kpep[i];
      sizi[threadIdx.x] = prepep[i];
      alphai[threadIdx.x] = dmppep[i];
      epli[threadIdx.x] = lpep[i];
      springk = kpep[k];
      sizk = prepep[k];
      alphak = dmppep[k];
      eplk = lpep[k];

      unsigned int dinfo0 = dinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dinfo0 & srcmask) == 0;
         real scaleb = 1;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if ((eplk or epli[klane]) and r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real ks2i[3][3], ks2k[3][3];
            pair_alterpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, springi[klane], sizi[klane],
               alphai[klane], springk, sizk, alphak, ks2i, ks2k);
            psci00[klane] += ks2i[0][0];
            psci01[klane] += ks2i[0][1];
            psci02[klane] += ks2i[0][2];
            psci10[klane] += ks2i[1][0];
            psci11[klane] += ks2i[1][1];
            psci12[klane] += ks2i[1][2];
            psci20[klane] += ks2i[2][0];
            psci21[klane] += ks2i[2][1];
            psci22[klane] += ks2i[2][2];
            psck00 += ks2k[0][0];
            psck01 += ks2k[0][1];
            psck02 += ks2k[0][2];
            psck10 += ks2k[1][0];
            psck11 += ks2k[1][1];
            psck12 += ks2k[1][2];
            psck20 += ks2k[2][0];
            psck21 += ks2k[2][1];
            psck22 += ks2k[2][2];
         }

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
      }

      atomic_add(psci00[threadIdx.x], &polscale[i][0]);
      atomic_add(psci01[threadIdx.x], &polscale[i][1]);
      atomic_add(psci02[threadIdx.x], &polscale[i][2]);
      atomic_add(psci10[threadIdx.x], &polscale[i][3]);
      atomic_add(psci11[threadIdx.x], &polscale[i][4]);
      atomic_add(psci12[threadIdx.x], &polscale[i][5]);
      atomic_add(psci20[threadIdx.x], &polscale[i][6]);
      atomic_add(psci21[threadIdx.x], &polscale[i][7]);
      atomic_add(psci22[threadIdx.x], &polscale[i][8]);
      atomic_add(psck00, &polscale[k][0]);
      atomic_add(psck01, &polscale[k][1]);
      atomic_add(psck02, &polscale[k][2]);
      atomic_add(psck10, &polscale[k][3]);
      atomic_add(psck11, &polscale[k][4]);
      atomic_add(psck12, &polscale[k][5]);
      atomic_add(psck20, &polscale[k][6]);
      atomic_add(psck21, &polscale[k][7]);
      atomic_add(psck22, &polscale[k][8]);
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      psci00[threadIdx.x] = 0;
      psci01[threadIdx.x] = 0;
      psci02[threadIdx.x] = 0;
      psci10[threadIdx.x] = 0;
      psci11[threadIdx.x] = 0;
      psci12[threadIdx.x] = 0;
      psci20[threadIdx.x] = 0;
      psci21[threadIdx.x] = 0;
      psci22[threadIdx.x] = 0;
      psck00 = 0;
      psck01 = 0;
      psck02 = 0;
      psck10 = 0;
      psck11 = 0;
      psck12 = 0;
      psck20 = 0;
      psck21 = 0;
      psck22 = 0;

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

      springi[threadIdx.x] = kpep[i];
      sizi[threadIdx.x] = prepep[i];
      alphai[threadIdx.x] = dmppep[i];
      epli[threadIdx.x] = lpep[i];
      springk = kpep[k];
      sizk = prepep[k];
      alphak = dmppep[k];
      eplk = lpep[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real scaleb = 1;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if ((eplk or epli[klane]) and r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real ks2i[3][3], ks2k[3][3];
            pair_alterpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, springi[klane], sizi[klane],
               alphai[klane], springk, sizk, alphak, ks2i, ks2k);
            psci00[klane] += ks2i[0][0];
            psci01[klane] += ks2i[0][1];
            psci02[klane] += ks2i[0][2];
            psci10[klane] += ks2i[1][0];
            psci11[klane] += ks2i[1][1];
            psci12[klane] += ks2i[1][2];
            psci20[klane] += ks2i[2][0];
            psci21[klane] += ks2i[2][1];
            psci22[klane] += ks2i[2][2];
            psck00 += ks2k[0][0];
            psck01 += ks2k[0][1];
            psck02 += ks2k[0][2];
            psck10 += ks2k[1][0];
            psck11 += ks2k[1][1];
            psck12 += ks2k[1][2];
            psck20 += ks2k[2][0];
            psck21 += ks2k[2][1];
            psck22 += ks2k[2][2];
         }
      }

      atomic_add(psci00[threadIdx.x], &polscale[i][0]);
      atomic_add(psci01[threadIdx.x], &polscale[i][1]);
      atomic_add(psci02[threadIdx.x], &polscale[i][2]);
      atomic_add(psci10[threadIdx.x], &polscale[i][3]);
      atomic_add(psci11[threadIdx.x], &polscale[i][4]);
      atomic_add(psci12[threadIdx.x], &polscale[i][5]);
      atomic_add(psci20[threadIdx.x], &polscale[i][6]);
      atomic_add(psci21[threadIdx.x], &polscale[i][7]);
      atomic_add(psci22[threadIdx.x], &polscale[i][8]);
      atomic_add(psck00, &polscale[k][0]);
      atomic_add(psck01, &polscale[k][1]);
      atomic_add(psck02, &polscale[k][2]);
      atomic_add(psck10, &polscale[k][3]);
      atomic_add(psck11, &polscale[k][4]);
      atomic_add(psck12, &polscale[k][5]);
      atomic_add(psck20, &polscale[k][6]);
      atomic_add(psck21, &polscale[k][7]);
      atomic_add(psck22, &polscale[k][8]);
   }
}

__global__
static void alterpolInit_cu1(int n, real (*restrict polscale)[3][3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      polscale[i][0][0] = 1.f;
      polscale[i][0][1] = 0.f;
      polscale[i][0][2] = 0.f;
      polscale[i][1][0] = 0.f;
      polscale[i][1][1] = 1.f;
      polscale[i][1][2] = 0.f;
      polscale[i][2][0] = 0.f;
      polscale[i][2][1] = 0.f;
      polscale[i][2][2] = 1.f;
   }
}

__global__
static void alterpolInvert_cu1(
   int n, real (*restrict polscale)[3][3], real (*restrict polinv)[3][3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real det;
      real(&ps)[3][3] = polscale[i];
      det = ps[0][0] * (ps[1][1] * ps[2][2] - ps[1][2] * ps[2][1]) -
         ps[1][0] * (ps[0][1] * ps[2][2] - ps[2][1] * ps[0][2]) +
         ps[2][0] * (ps[0][1] * ps[1][2] - ps[1][1] * ps[0][2]);
      polinv[i][0][0] = (ps[1][1] * ps[2][2] - ps[1][2] * ps[2][1]) / det;
      polinv[i][1][0] = (ps[2][0] * ps[1][2] - ps[1][0] * ps[2][2]) / det;
      polinv[i][2][0] = (ps[1][0] * ps[2][1] - ps[2][0] * ps[1][1]) / det;
      polinv[i][0][1] = (ps[2][1] * ps[0][2] - ps[0][1] * ps[2][2]) / det;
      polinv[i][1][1] = (ps[0][0] * ps[2][2] - ps[2][0] * ps[0][2]) / det;
      polinv[i][2][1] = (ps[0][1] * ps[2][0] - ps[0][0] * ps[2][1]) / det;
      polinv[i][0][2] = (ps[0][1] * ps[1][2] - ps[0][2] * ps[1][1]) / det;
      polinv[i][1][2] = (ps[0][2] * ps[1][0] - ps[0][0] * ps[1][2]) / det;
      polinv[i][2][2] = (ps[0][0] * ps[1][1] - ps[0][1] * ps[1][0]) / det;
   }
}

void alterpol_cu(real (*polscale)[3][3], real (*polinv)[3][3])
{
   const auto& st = *mspatial_v2_unit;
   real cut = switchCut(Switch::REPULS);
   real off = switchOff(Switch::REPULS);

   launch_k1s(g::s0, n, alterpolInit_cu1, //
      n, polscale);

   int ngrid = gpuGridSize(BLOCK_DIM);
   alterpol_cu1<<<ngrid, BLOCK_DIM, 0, g::s0>>>(n, TINKER_IMAGE_ARGS, cut, off, st.si2.bit0,
      nmdwexclude, mdwexclude, mdwexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
      st.niak, st.iak, st.lst, reinterpret_cast<real(*)[9]>(polscale), kpep, prepep, dmppep, lpep,
      scrtyp);

   launch_k1s(g::s0, n, alterpolInvert_cu1, //
      n, polscale, polinv);
}

// ck.py Version 2.0.3

template <class Ver>
__global__
void dexpol_cu1(int n, TINKER_IMAGE_PARAMS, VirialBuffer restrict vep, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, real cut, real off,
   const unsigned* restrict dinfo, int nexclude, const int (*restrict exclude)[2],
   const real (*restrict exclude_scale)[3], const real* restrict x, const real* restrict y,
   const real* restrict z, const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   const real* restrict polarity, const real (*restrict uind)[3], const real* restrict kpep,
   const real* restrict prepep, const real* restrict dmppep, const int* restrict lpep,
   ExpolScr scrtyp, real f)
{
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec veptlxx, veptlyx, veptlzx, veptlyy, veptlzy, veptlzz;
   if CONSTEXPR (do_v) {
      veptlxx = 0;
      veptlyx = 0;
      veptlzx = 0;
      veptlyy = 0;
      veptlzy = 0;
      veptlzz = 0;
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
   real frcxk;
   real frcyk;
   real frczk;
   __shared__ real uix[BLOCK_DIM];
   __shared__ real uiy[BLOCK_DIM];
   __shared__ real uiz[BLOCK_DIM];
   __shared__ real springi[BLOCK_DIM];
   __shared__ real sizi[BLOCK_DIM];
   __shared__ real alphai[BLOCK_DIM];
   __shared__ int epli[BLOCK_DIM];
   __shared__ real poli[BLOCK_DIM];
   real ukx;
   real uky;
   real ukz;
   real springk;
   real sizk;
   real alphak;
   int eplk;
   real polk;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      frcxi[threadIdx.x] = 0;
      frcyi[threadIdx.x] = 0;
      frczi[threadIdx.x] = 0;
      frcxk = 0;
      frcyk = 0;
      frczk = 0;

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];

      xi[klane] = x[i];
      yi[klane] = y[i];
      zi[klane] = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      uix[klane] = uind[i][0];
      uiy[klane] = uind[i][1];
      uiz[klane] = uind[i][2];
      springi[klane] = kpep[i];
      sizi[klane] = prepep[i];
      alphai[klane] = dmppep[i];
      epli[klane] = lpep[i];
      poli[klane] = polarity[i];
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      springk = kpep[k];
      sizk = prepep[k];
      alphak = dmppep[k];
      eplk = lpep[k];
      polk = polarity[k];

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];
      real r2 = image2(xr, yr, zr);
      if ((eplk or epli[klane]) and r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real frc[3];
         pair_dexpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, uix[klane], uiy[klane], uiz[klane],
            ukx, uky, ukz, springi[klane] / poli[klane], sizi[klane], alphai[klane], springk / polk,
            sizk, alphak, f, frc);
         frcxi[klane] += frc[0];
         frcyi[klane] += frc[1];
         frczi[klane] += frc[2];
         frcxk -= frc[0];
         frcyk -= frc[1];
         frczk -= frc[2];

         if CONSTEXPR (do_v) {
            real vxx = -xr * frc[0];
            real vxy = -0.5f * (yr * frc[0] + xr * frc[1]);
            real vxz = -0.5f * (zr * frc[0] + xr * frc[2]);
            real vyy = -yr * frc[1];
            real vyz = -0.5f * (zr * frc[1] + yr * frc[2]);
            real vzz = -zr * frc[2];
            veptlxx += floatTo<vbuf_prec>(vxx);
            veptlyx += floatTo<vbuf_prec>(vxy);
            veptlzx += floatTo<vbuf_prec>(vxz);
            veptlyy += floatTo<vbuf_prec>(vyy);
            veptlzy += floatTo<vbuf_prec>(vyz);
            veptlzz += floatTo<vbuf_prec>(vzz);
         }
      }

      atomic_add(frcxi[threadIdx.x], gx, i);
      atomic_add(frcyi[threadIdx.x], gy, i);
      atomic_add(frczi[threadIdx.x], gz, i);
      atomic_add(frcxk, gx, k);
      atomic_add(frcyk, gy, k);
      atomic_add(frczk, gz, k);
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      frcxi[threadIdx.x] = 0;
      frcyi[threadIdx.x] = 0;
      frczi[threadIdx.x] = 0;
      frcxk = 0;
      frcyk = 0;
      frczk = 0;

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

      uix[threadIdx.x] = uind[i][0];
      uiy[threadIdx.x] = uind[i][1];
      uiz[threadIdx.x] = uind[i][2];
      springi[threadIdx.x] = kpep[i];
      sizi[threadIdx.x] = prepep[i];
      alphai[threadIdx.x] = dmppep[i];
      epli[threadIdx.x] = lpep[i];
      poli[threadIdx.x] = polarity[i];
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      springk = kpep[k];
      sizk = prepep[k];
      alphak = dmppep[k];
      eplk = lpep[k];
      polk = polarity[k];

      unsigned int dinfo0 = dinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dinfo0 & srcmask) == 0;
         real scaleb = 1;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if ((eplk or epli[klane]) and r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real frc[3];
            pair_dexpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, uix[klane], uiy[klane], uiz[klane],
               ukx, uky, ukz, springi[klane] / poli[klane], sizi[klane], alphai[klane],
               springk / polk, sizk, alphak, f, frc);
            frcxi[klane] += frc[0];
            frcyi[klane] += frc[1];
            frczi[klane] += frc[2];
            frcxk -= frc[0];
            frcyk -= frc[1];
            frczk -= frc[2];

            if CONSTEXPR (do_v) {
               real vxx = -xr * frc[0];
               real vxy = -0.5f * (yr * frc[0] + xr * frc[1]);
               real vxz = -0.5f * (zr * frc[0] + xr * frc[2]);
               real vyy = -yr * frc[1];
               real vyz = -0.5f * (zr * frc[1] + yr * frc[2]);
               real vzz = -zr * frc[2];
               veptlxx += floatTo<vbuf_prec>(vxx);
               veptlyx += floatTo<vbuf_prec>(vxy);
               veptlzx += floatTo<vbuf_prec>(vxz);
               veptlyy += floatTo<vbuf_prec>(vyy);
               veptlzy += floatTo<vbuf_prec>(vyz);
               veptlzz += floatTo<vbuf_prec>(vzz);
            }
         }

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
      }

      atomic_add(frcxi[threadIdx.x], gx, i);
      atomic_add(frcyi[threadIdx.x], gy, i);
      atomic_add(frczi[threadIdx.x], gz, i);
      atomic_add(frcxk, gx, k);
      atomic_add(frcyk, gy, k);
      atomic_add(frczk, gz, k);
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      frcxi[threadIdx.x] = 0;
      frcyi[threadIdx.x] = 0;
      frczi[threadIdx.x] = 0;
      frcxk = 0;
      frcyk = 0;
      frczk = 0;

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

      uix[threadIdx.x] = uind[i][0];
      uiy[threadIdx.x] = uind[i][1];
      uiz[threadIdx.x] = uind[i][2];
      springi[threadIdx.x] = kpep[i];
      sizi[threadIdx.x] = prepep[i];
      alphai[threadIdx.x] = dmppep[i];
      epli[threadIdx.x] = lpep[i];
      poli[threadIdx.x] = polarity[i];
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      springk = kpep[k];
      sizk = prepep[k];
      alphak = dmppep[k];
      eplk = lpep[k];
      polk = polarity[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real scaleb = 1;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if ((eplk or epli[klane]) and r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real frc[3];
            pair_dexpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, uix[klane], uiy[klane], uiz[klane],
               ukx, uky, ukz, springi[klane] / poli[klane], sizi[klane], alphai[klane],
               springk / polk, sizk, alphak, f, frc);
            frcxi[klane] += frc[0];
            frcyi[klane] += frc[1];
            frczi[klane] += frc[2];
            frcxk -= frc[0];
            frcyk -= frc[1];
            frczk -= frc[2];

            if CONSTEXPR (do_v) {
               real vxx = -xr * frc[0];
               real vxy = -0.5f * (yr * frc[0] + xr * frc[1]);
               real vxz = -0.5f * (zr * frc[0] + xr * frc[2]);
               real vyy = -yr * frc[1];
               real vyz = -0.5f * (zr * frc[1] + yr * frc[2]);
               real vzz = -zr * frc[2];
               veptlxx += floatTo<vbuf_prec>(vxx);
               veptlyx += floatTo<vbuf_prec>(vxy);
               veptlzx += floatTo<vbuf_prec>(vxz);
               veptlyy += floatTo<vbuf_prec>(vyy);
               veptlzy += floatTo<vbuf_prec>(vyz);
               veptlzz += floatTo<vbuf_prec>(vzz);
            }
         }
      }

      atomic_add(frcxi[threadIdx.x], gx, i);
      atomic_add(frcyi[threadIdx.x], gy, i);
      atomic_add(frczi[threadIdx.x], gz, i);
      atomic_add(frcxk, gx, k);
      atomic_add(frcyk, gy, k);
      atomic_add(frczk, gz, k);
   }

   if CONSTEXPR (do_v) {
      atomic_add(veptlxx, veptlyx, veptlzx, veptlyy, veptlzy, veptlzz, vep, ithread);
   }
}

void dexpol_cu(int vers, const real (*uind)[3], grad_prec* depx, grad_prec* depy, grad_prec* depz,
   VirialBuffer restrict vir_ep)
{
   const auto& st = *mspatial_v2_unit;
   real cut = switchCut(Switch::REPULS);
   real off = switchOff(Switch::REPULS);

   real f = 0.5f * electric / dielec;

   int ngrid = gpuGridSize(BLOCK_DIM);

#define DEXPOL_CU1_ARGS                                                                            \
   n, TINKER_IMAGE_ARGS, vir_ep, depx, depy, depz, cut, off, st.si2.bit0, nmdwexclude, mdwexclude, \
      mdwexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst,  \
      polarity, uind, kpep, prepep, dmppep, lpep, scrtyp, f

   if (vers & calc::virial) {
      dexpol_cu1<calc::V6><<<ngrid, BLOCK_DIM, 0, g::s0>>>(DEXPOL_CU1_ARGS);
   } else if (vers & calc::grad) {
      dexpol_cu1<calc::V5><<<ngrid, BLOCK_DIM, 0, g::s0>>>(DEXPOL_CU1_ARGS);
   } else {
      assert(false &&
         "This function should not have been called if neither gradient nor virial is calculated.");
   }

#undef DEXPOL_CU1_ARGS
}
}

#include "ff/amoeba/induce.h"
#include "ff/amoebamod.h"
#include "ff/atom.h"
#include "ff/hippo/induce.h"
#include "ff/hippomod.h"
#include "ff/switch.h"
#include "seq/launch.h"
#include "tool/darray.h"
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>

namespace tinker {
__global__
void eppcgUdirDonly(
   int n, const real* restrict polarity, real (*restrict udir)[3], const real (*restrict field)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli = polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         udir[i][j] = poli * field[i][j];
      }
   }
}

__global__
void eppcgUdirGuess(int n, const real* restrict polarity, real (*restrict uind)[3],
   const real (*restrict field)[3], const real (*restrict polinv)[3][3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli = polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         uind[i][j] = poli *
            (polinv[i][0][j] * field[i][0] + polinv[i][1][j] * field[i][1] +
               polinv[i][2][j] * field[i][2]);
      }
   }
}

__global__
void eppcgRsd2(int n, const real* restrict polarity_inv, //
   real (*restrict rsd)[3],                              //
   const real (*restrict udir)[3], const real (*restrict uind)[3], const real (*restrict field)[3],
   const real (*restrict polscale)[3][3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         rsd[i][j] = (udir[i][j] - uind[i][0] * polscale[i][0][j] - uind[i][1] * polscale[i][1][j] -
                        uind[i][2] * polscale[i][2][j]) *
               poli_inv +
            field[i][j];
      }
   }
}

__global__
void eppcgRsd1(int n, const real* restrict polarity, real (*restrict rsd)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      if (polarity[i] == 0) {
         rsd[i][0] = 0;
         rsd[i][1] = 0;
         rsd[i][2] = 0;
      }
   }
}

__global__
void eppcgP4(int n, const real* restrict polarity_inv, real (*restrict vec)[3],
   const real (*restrict conj)[3], const real (*restrict field)[3],
   const real (*restrict polscale)[3][3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         vec[i][j] = poli_inv *
               (conj[i][0] * polscale[i][0][j] + conj[i][1] * polscale[i][1][j] +
                  conj[i][2] * polscale[i][2][j]) -
            field[i][j];
   }
}

__global__
void eppcgP5(int n, const real* restrict polarity, //
   const real* restrict ka,                        //
   const real* restrict ksum, real (*restrict uind)[3], const real (*restrict conj)[3],
   real (*restrict rsd)[3], const real (*restrict vec)[3])
{
   real kaval = *ka;
   real a = *ksum / kaval;
   if (kaval == 0) a = 0;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         uind[i][j] += a * conj[i][j];
         rsd[i][j] -= a * vec[i][j];
      }
      if (polarity[i] == 0) {
         rsd[i][0] = 0;
         rsd[i][1] = 0;
         rsd[i][2] = 0;
      }
   }
}

__global__
void eppcgP6(int n, const real* restrict ksum, const real* restrict ksum1, real (*restrict conj)[3],
   real (*restrict zrsd)[3])
{
   real ksumval = *ksum;
   real b = *ksum1 / ksumval;
   if (ksumval == 0) b = 0;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         conj[i][j] = zrsd[i][j] + b * conj[i][j];
   }
}

__global__
void eppcgPeek1(int n, float pcgpeek, const real* restrict polarity, real (*restrict uind)[3],
   const real (*restrict rsd)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real term = pcgpeek * polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         uind[i][j] += term * rsd[i][j];
   }
}

void induceMutualPcg4_cu(real (*uind)[3])
{
   auto* field = work01_;
   auto* rsd = work02_;
   auto* zrsd = work03_;
   auto* conj = work04_;
   auto* vec = work05_;

   const bool sparse_prec = polpcg::pcgprec and (switchOff(Switch::USOLVE) > 0);
   bool dirguess = polpcg::pcgguess;
   bool predict = polpred != UPred::NONE;
   if (predict and nualt < maxualt) {
      predict = false;
      dirguess = true;
   }

   // get the electrostatic field due to permanent multipoles
   dfieldChgpen(field);
   // direct induced dipoles
   launch_k1s(g::s0, n, eppcgUdirDonly, n, polarity, udir, field);

   alterpol(polscale, polinv);

   // initial induced dipole
   if (predict) {
      ulspredSum(uind, nullptr);
   } else if (dirguess) {
      launch_k1s(g::s0, n, eppcgUdirGuess, n, polarity, uind, field, polinv);
   } else {
      darray::zero(g::q0, n, uind);
   }

   if (predict) {
      ufieldChgpen(uind, field);
      launch_k1s(g::s0, n, eppcgRsd2, n, polarity_inv, rsd, udir, uind, field, polscale);
   } else if (dirguess) {
      // uind is used here instead of udir since without exchange polarization udir = uind
      // but with exchange polarization udir != uind (for dirguess).
      ufieldChgpen(uind, rsd);
   } else {
      darray::copy(g::q0, n, rsd, field);
   }
   launch_k1s(g::s0, n, eppcgRsd1, n, polarity, rsd);

   // initial M r(0) and p(0)
   if (sparse_prec) {
      sparsePrecondBuild2();
      sparsePrecondApply2(rsd, zrsd);
   } else {
      diagPrecond2(rsd, zrsd);
   }
   darray::copy(g::q0, n, conj, zrsd);

   // initial r(0) M r(0)
   real* sum = &((real*)dptr_buf)[0];
   darray::dot(g::q0, n, sum, rsd, zrsd);

   // conjugate gradient iteration of the mutual induced dipoles
   const bool debug = inform::debug;
   const int politer = polpot::politer;
   const real poleps = polpot::poleps;
   const real debye = units::debye;
   const real pcgpeek = polpcg::pcgpeek;
   const int maxiter = 100; // see also subroutine induce0a in induce.f
   const int miniter = std::min(3, n);

   bool done = false;
   int iter = 0;
   real eps = 100;
   // real epsold;

   while (not done) {
      ++iter;

      // T p and p
      // vec = (inv_alpha + Tu) conj, field = -Tu conj
      // vec = inv_alpha * conj - field
      ufieldChgpen(conj, field);
      launch_k1s(g::s0, n, eppcgP4, n, polarity_inv, vec, conj, field, polscale);

      // a <- p T p
      real* a = &((real*)dptr_buf)[1];
      // a <- r M r / p T p; a = sum / a; ap = sump / ap
      darray::dot(g::q0, n, a, conj, vec);

      // u <- u + a p
      // r <- r - a T p
      launch_k1s(g::s0, n, eppcgP5, n, polarity, a, sum, uind, conj, rsd, vec);

      // calculate/update M r
      if (sparse_prec)
         sparsePrecondApply2(rsd, zrsd);
      else
         diagPrecond2(rsd, zrsd);

      // b = sum1 / sum; bp = sump1 / sump
      real* sum1 = &((real*)dptr_buf)[2];
      darray::dot(g::q0, n, sum1, rsd, zrsd);

      // calculate/update p
      launch_k1s(g::s0, n, eppcgP6, n, sum, sum1, conj, zrsd);

      // copy sum1/p to sum/p
      darray::copy(g::q0, 2, sum, sum1);

      real* epsd = &((real*)dptr_buf)[3];
      darray::dot(g::q0, n, epsd, rsd, rsd);
      check_rt(
         cudaMemcpyAsync((real*)pinned_buf, epsd, sizeof(real), cudaMemcpyDeviceToHost, g::s0));
      check_rt(cudaStreamSynchronize(g::s0));
      // epsold = eps;
      eps = ((real*)pinned_buf)[0];
      eps = debye * REAL_SQRT(eps / n);

      if (debug) {
         if (iter == 1) {
            print(stdout,
               "\n Determination of SCF Induced Dipole Moments\n\n"
               "    Iter    RMS Residual (Debye)\n\n");
         }
         print(stdout, " %8d       %-16.10f\n", iter, eps);
      }

      if (eps < poleps) done = true;
      // if (eps > epsold) done = true;
      if (iter < miniter) done = false;
      if (iter >= politer) done = true;

      // apply a "peek" iteration to the mutual induced dipoles
      if (done) launch_k1s(g::s0, n, eppcgPeek1, n, pcgpeek, polarity, uind, rsd);
   }

   // print the results from the conjugate gradient iteration
   if (debug) {
      print(stdout,
         " Induced Dipoles :    Iterations %4d      RMS"
         " Residual %14.10f\n",
         iter, eps);
   }

   // terminate the calculation if dipoles failed to converge
   if (iter >= maxiter) {
      printError();
      TINKER_THROW("INDUCE  --  Warning, Induced Dipoles are not Converged");
   }
}
}
