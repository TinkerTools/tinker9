// ck.py Version 3.0.2
__global__
void sparsePrecond_cu4(int n, TINKER_IMAGE_PARAMS, real off, const unsigned* restrict winfo, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale, const real* restrict x, const real* restrict y,
   const real* restrict z, const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, const real (*restrict rsd)[3], real (*restrict zrsd)[3],
   const real* restrict palpha, const real* restrict polarity)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   __shared__ real uidx[BLOCK_DIM], uidy[BLOCK_DIM], uidz[BLOCK_DIM], alphai[BLOCK_DIM], poli[BLOCK_DIM];
   real xi, yi, zi;
   real xk, yk, zk, ukdx, ukdy, ukdz, alphak, polk;
   real fidx, fidy, fidz;
   real fkdx, fkdy, fkdz;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      fidx = 0;
      fidy = 0;
      fidz = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii];

      uidx[klane] = rsd[i][0];
      uidy[klane] = rsd[i][1];
      uidz[klane] = rsd[i][2];
      alphai[klane] = palpha[i];
      poli[klane] = polarity[i];
      xi = x[i];
      yi = y[i];
      zi = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      alphak = palpha[k];
      polk = polarity[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real dmpik[3];
         damp_mut(dmpik, r, alphai[klane], alphak);
         real scale3 = scalea * dmpik[1];
         real scale5 = scalea * dmpik[2];

         real polik = poli[klane] * polk;
         real rr3 = scale3 * polik * REAL_RECIP(r * r2);
         real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

         real c;
         c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
         fidx += c * xr - rr3 * ukdx;
         fidy += c * yr - rr3 * ukdy;
         fidz += c * zr - rr3 * ukdz;

         c = rr5 * dot3(xr, yr, zr, uidx[klane], uidy[klane], uidz[klane]);
         fkdx += c * xr - rr3 * uidx[klane];
         fkdy += c * yr - rr3 * uidy[klane];
         fkdz += c * zr - rr3 * uidz[klane];
      } // end if (include)

      atomic_add(fidx, &zrsd[i][0]);
      atomic_add(fidy, &zrsd[i][1]);
      atomic_add(fidz, &zrsd[i][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      fidx = 0;
      fidy = 0;
      fidz = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;

      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);

      int iid = ty * WARP_SIZE + ilane;
      int atomi = min(iid, n - 1);
      int i = sorted[atomi].unsorted;
      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      int k = sorted[atomk].unsorted;
      uidx[threadIdx.x] = rsd[i][0];
      uidy[threadIdx.x] = rsd[i][1];
      uidz[threadIdx.x] = rsd[i][2];
      alphai[threadIdx.x] = palpha[i];
      poli[threadIdx.x] = polarity[i];
      xi = sorted[atomi].x;
      yi = sorted[atomi].y;
      zi = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      alphak = palpha[k];
      polk = polarity[k];
      __syncwarp();

      unsigned int winfo0 = winfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (winfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real dmpik[3];
            damp_mut(dmpik, r, alphai[klane], alphak);
            real scale3 = scalea * dmpik[1];
            real scale5 = scalea * dmpik[2];

            real polik = poli[klane] * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

            real c;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            fidx += c * xr - rr3 * ukdx;
            fidy += c * yr - rr3 * ukdy;
            fidz += c * zr - rr3 * ukdz;

            c = rr5 * dot3(xr, yr, zr, uidx[klane], uidy[klane], uidz[klane]);
            fkdx += c * xr - rr3 * uidx[klane];
            fkdy += c * yr - rr3 * uidy[klane];
            fkdz += c * zr - rr3 * uidz[klane];
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         fidx = __shfl_sync(ALL_LANES, fidx, ilane + 1);
         fidy = __shfl_sync(ALL_LANES, fidy, ilane + 1);
         fidz = __shfl_sync(ALL_LANES, fidz, ilane + 1);
      }

      atomic_add(fidx, &zrsd[i][0]);
      atomic_add(fidy, &zrsd[i][1]);
      atomic_add(fidz, &zrsd[i][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
      __syncwarp();
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      fidx = 0;
      fidy = 0;
      fidz = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;

      int ty = iak[iw];
      int atomi = ty * WARP_SIZE + ilane;
      int i = sorted[atomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      uidx[threadIdx.x] = rsd[i][0];
      uidy[threadIdx.x] = rsd[i][1];
      uidz[threadIdx.x] = rsd[i][2];
      alphai[threadIdx.x] = palpha[i];
      poli[threadIdx.x] = polarity[i];
      xi = sorted[atomi].x;
      yi = sorted[atomi].y;
      zi = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      alphak = palpha[k];
      polk = polarity[k];
      __syncwarp();

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real dmpik[3];
            damp_mut(dmpik, r, alphai[klane], alphak);
            real scale3 = scalea * dmpik[1];
            real scale5 = scalea * dmpik[2];

            real polik = poli[klane] * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

            real c;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            fidx += c * xr - rr3 * ukdx;
            fidy += c * yr - rr3 * ukdy;
            fidz += c * zr - rr3 * ukdz;

            c = rr5 * dot3(xr, yr, zr, uidx[klane], uidy[klane], uidz[klane]);
            fkdx += c * xr - rr3 * uidx[klane];
            fkdy += c * yr - rr3 * uidy[klane];
            fkdz += c * zr - rr3 * uidz[klane];
         } // end if (include)

         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         fidx = __shfl_sync(ALL_LANES, fidx, ilane + 1);
         fidy = __shfl_sync(ALL_LANES, fidy, ilane + 1);
         fidz = __shfl_sync(ALL_LANES, fidz, ilane + 1);
      }

      atomic_add(fidx, &zrsd[i][0]);
      atomic_add(fidy, &zrsd[i][1]);
      atomic_add(fidz, &zrsd[i][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
      __syncwarp();
   }
}
