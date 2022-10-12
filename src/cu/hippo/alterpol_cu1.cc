// ck.py Version 3.0.2
__global__
void alterpol_cu1(int n, TINKER_IMAGE_PARAMS, real cut, real off, const unsigned* restrict dinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3], const real* restrict x,
   const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst, real (*restrict polscale)[9],
   const real* restrict kpep, const real* restrict prepep, const real* restrict dmppep, const int* restrict lpep,
   ExpolScr scrtyp)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   __shared__ real xi[BLOCK_DIM], yi[BLOCK_DIM], zi[BLOCK_DIM], springi[BLOCK_DIM], sizi[BLOCK_DIM], alphai[BLOCK_DIM];
   __shared__ int epli[BLOCK_DIM];
   __shared__ real xk[BLOCK_DIM], yk[BLOCK_DIM], zk[BLOCK_DIM], springk[BLOCK_DIM], sizk[BLOCK_DIM], alphak[BLOCK_DIM];
   __shared__ int eplk[BLOCK_DIM];
   real psci00, psci01, psci02, psci10, psci11, psci12, psci20, psci21, psci22;
   real psck00, psck01, psck02, psck10, psck11, psck12, psck20, psck21, psck22;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      psci00 = 0;
      psci01 = 0;
      psci02 = 0;
      psci10 = 0;
      psci11 = 0;
      psci12 = 0;
      psci20 = 0;
      psci21 = 0;
      psci22 = 0;
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
      springi[klane] = kpep[i];
      sizi[klane] = prepep[i];
      alphai[klane] = dmppep[i];
      epli[klane] = lpep[i];
      xk[threadIdx.x] = x[k];
      yk[threadIdx.x] = y[k];
      zk[threadIdx.x] = z[k];
      springk[threadIdx.x] = kpep[k];
      sizk[threadIdx.x] = prepep[k];
      alphak[threadIdx.x] = dmppep[k];
      eplk[threadIdx.x] = lpep[k];

      constexpr bool incl = true;
      real xr = xk[threadIdx.x] - xi[klane];
      real yr = yk[threadIdx.x] - yi[klane];
      real zr = zk[threadIdx.x] - zi[klane];
      real r2 = image2(xr, yr, zr);
      if ((eplk[threadIdx.x] or epli[klane]) and r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real ks2i[3][3], ks2k[3][3];
         pair_alterpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, springi[klane], sizi[klane], alphai[klane],
            springk[threadIdx.x], sizk[threadIdx.x], alphak[threadIdx.x], ks2i, ks2k);
         psci00 += ks2i[0][0];
         psci01 += ks2i[0][1];
         psci02 += ks2i[0][2];
         psci10 += ks2i[1][0];
         psci11 += ks2i[1][1];
         psci12 += ks2i[1][2];
         psci20 += ks2i[2][0];
         psci21 += ks2i[2][1];
         psci22 += ks2i[2][2];
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

      atomic_add(psci00, &polscale[i][0]);
      atomic_add(psci01, &polscale[i][1]);
      atomic_add(psci02, &polscale[i][2]);
      atomic_add(psci10, &polscale[i][3]);
      atomic_add(psci11, &polscale[i][4]);
      atomic_add(psci12, &polscale[i][5]);
      atomic_add(psci20, &polscale[i][6]);
      atomic_add(psci21, &polscale[i][7]);
      atomic_add(psci22, &polscale[i][8]);
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
      psci00 = 0;
      psci01 = 0;
      psci02 = 0;
      psci10 = 0;
      psci11 = 0;
      psci12 = 0;
      psci20 = 0;
      psci21 = 0;
      psci22 = 0;
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
      springi[threadIdx.x] = kpep[i];
      sizi[threadIdx.x] = prepep[i];
      alphai[threadIdx.x] = dmppep[i];
      epli[threadIdx.x] = lpep[i];
      xk[threadIdx.x] = sorted[atomk].x;
      yk[threadIdx.x] = sorted[atomk].y;
      zk[threadIdx.x] = sorted[atomk].z;
      springk[threadIdx.x] = kpep[k];
      sizk[threadIdx.x] = prepep[k];
      alphak[threadIdx.x] = dmppep[k];
      eplk[threadIdx.x] = lpep[k];
      __syncwarp();

      unsigned int dinfo0 = dinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dinfo0 & srcmask) == 0;
         real scaleb = 1;
         real xr = xk[threadIdx.x] - xi[klane];
         real yr = yk[threadIdx.x] - yi[klane];
         real zr = zk[threadIdx.x] - zi[klane];
         real r2 = image2(xr, yr, zr);
         if ((eplk[threadIdx.x] or epli[klane]) and r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real ks2i[3][3], ks2k[3][3];
            pair_alterpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, springi[klane], sizi[klane], alphai[klane],
               springk[threadIdx.x], sizk[threadIdx.x], alphak[threadIdx.x], ks2i, ks2k);
            psci00 += ks2i[0][0];
            psci01 += ks2i[0][1];
            psci02 += ks2i[0][2];
            psci10 += ks2i[1][0];
            psci11 += ks2i[1][1];
            psci12 += ks2i[1][2];
            psci20 += ks2i[2][0];
            psci21 += ks2i[2][1];
            psci22 += ks2i[2][2];
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
         psci00 = __shfl_sync(ALL_LANES, psci00, ilane + 1);
         psci01 = __shfl_sync(ALL_LANES, psci01, ilane + 1);
         psci02 = __shfl_sync(ALL_LANES, psci02, ilane + 1);
         psci10 = __shfl_sync(ALL_LANES, psci10, ilane + 1);
         psci11 = __shfl_sync(ALL_LANES, psci11, ilane + 1);
         psci12 = __shfl_sync(ALL_LANES, psci12, ilane + 1);
         psci20 = __shfl_sync(ALL_LANES, psci20, ilane + 1);
         psci21 = __shfl_sync(ALL_LANES, psci21, ilane + 1);
         psci22 = __shfl_sync(ALL_LANES, psci22, ilane + 1);
      }

      atomic_add(psci00, &polscale[i][0]);
      atomic_add(psci01, &polscale[i][1]);
      atomic_add(psci02, &polscale[i][2]);
      atomic_add(psci10, &polscale[i][3]);
      atomic_add(psci11, &polscale[i][4]);
      atomic_add(psci12, &polscale[i][5]);
      atomic_add(psci20, &polscale[i][6]);
      atomic_add(psci21, &polscale[i][7]);
      atomic_add(psci22, &polscale[i][8]);
      atomic_add(psck00, &polscale[k][0]);
      atomic_add(psck01, &polscale[k][1]);
      atomic_add(psck02, &polscale[k][2]);
      atomic_add(psck10, &polscale[k][3]);
      atomic_add(psck11, &polscale[k][4]);
      atomic_add(psck12, &polscale[k][5]);
      atomic_add(psck20, &polscale[k][6]);
      atomic_add(psck21, &polscale[k][7]);
      atomic_add(psck22, &polscale[k][8]);
      __syncwarp();
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      psci00 = 0;
      psci01 = 0;
      psci02 = 0;
      psci10 = 0;
      psci11 = 0;
      psci12 = 0;
      psci20 = 0;
      psci21 = 0;
      psci22 = 0;
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
      springi[threadIdx.x] = kpep[i];
      sizi[threadIdx.x] = prepep[i];
      alphai[threadIdx.x] = dmppep[i];
      epli[threadIdx.x] = lpep[i];
      xk[threadIdx.x] = sorted[atomk].x;
      yk[threadIdx.x] = sorted[atomk].y;
      zk[threadIdx.x] = sorted[atomk].z;
      springk[threadIdx.x] = kpep[k];
      sizk[threadIdx.x] = prepep[k];
      alphak[threadIdx.x] = dmppep[k];
      eplk[threadIdx.x] = lpep[k];
      __syncwarp();

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real scaleb = 1;
         real xr = xk[threadIdx.x] - xi[klane];
         real yr = yk[threadIdx.x] - yi[klane];
         real zr = zk[threadIdx.x] - zi[klane];
         real r2 = image2(xr, yr, zr);
         if ((eplk[threadIdx.x] or epli[klane]) and r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real ks2i[3][3], ks2k[3][3];
            pair_alterpol(scrtyp, r, scaleb, cut, off, xr, yr, zr, springi[klane], sizi[klane], alphai[klane],
               springk[threadIdx.x], sizk[threadIdx.x], alphak[threadIdx.x], ks2i, ks2k);
            psci00 += ks2i[0][0];
            psci01 += ks2i[0][1];
            psci02 += ks2i[0][2];
            psci10 += ks2i[1][0];
            psci11 += ks2i[1][1];
            psci12 += ks2i[1][2];
            psci20 += ks2i[2][0];
            psci21 += ks2i[2][1];
            psci22 += ks2i[2][2];
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

         psci00 = __shfl_sync(ALL_LANES, psci00, ilane + 1);
         psci01 = __shfl_sync(ALL_LANES, psci01, ilane + 1);
         psci02 = __shfl_sync(ALL_LANES, psci02, ilane + 1);
         psci10 = __shfl_sync(ALL_LANES, psci10, ilane + 1);
         psci11 = __shfl_sync(ALL_LANES, psci11, ilane + 1);
         psci12 = __shfl_sync(ALL_LANES, psci12, ilane + 1);
         psci20 = __shfl_sync(ALL_LANES, psci20, ilane + 1);
         psci21 = __shfl_sync(ALL_LANES, psci21, ilane + 1);
         psci22 = __shfl_sync(ALL_LANES, psci22, ilane + 1);
      }

      atomic_add(psci00, &polscale[i][0]);
      atomic_add(psci01, &polscale[i][1]);
      atomic_add(psci02, &polscale[i][2]);
      atomic_add(psci10, &polscale[i][3]);
      atomic_add(psci11, &polscale[i][4]);
      atomic_add(psci12, &polscale[i][5]);
      atomic_add(psci20, &polscale[i][6]);
      atomic_add(psci21, &polscale[i][7]);
      atomic_add(psci22, &polscale[i][8]);
      atomic_add(psck00, &polscale[k][0]);
      atomic_add(psck01, &polscale[k][1]);
      atomic_add(psck02, &polscale[k][2]);
      atomic_add(psck10, &polscale[k][3]);
      atomic_add(psck11, &polscale[k][4]);
      atomic_add(psck12, &polscale[k][5]);
      atomic_add(psck20, &polscale[k][6]);
      atomic_add(psck21, &polscale[k][7]);
      atomic_add(psck22, &polscale[k][8]);
      __syncwarp();
   }
}
