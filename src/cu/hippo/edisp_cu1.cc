// ck.py Version 3.1.0
template <class Ver, class DTYP>
__global__
void edisp_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nd, EnergyBuffer restrict ed, VirialBuffer restrict vd,
   grad_prec* restrict gx, grad_prec* restrict gy, grad_prec* restrict gz, real cut, real off,
   const unsigned* restrict dinfo, int nexclude, const int (*restrict exclude)[2], const real* restrict exclude_scale,
   const real* restrict x, const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted,
   int nakpl, const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   const real* restrict csix, const real* restrict adisp, real aewald, const int* restrict mut, real vlam, Vdw vcouple)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int ndtl;
   if CONSTEXPR (do_a) {
      ndtl = 0;
   }
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec edtl;
   if CONSTEXPR (do_e) {
      edtl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec vdtlxx, vdtlyx, vdtlzx, vdtlyy, vdtlzy, vdtlzz;
   if CONSTEXPR (do_v) {
      vdtlxx = 0;
      vdtlyx = 0;
      vdtlzx = 0;
      vdtlyy = 0;
      vdtlzy = 0;
      vdtlzz = 0;
   }
   real xi, yi, zi, ci, ai;
   int imut;
   real xk, yk, zk, ck, ak;
   int kmut;
   real gxi, gyi, gzi;
   real gxk, gyk, gzk;

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
      real scalea = exclude_scale[ii];

      xi = x[i];
      yi = y[i];
      zi = z[i];
      ci = csix[i];
      ai = adisp[i];
      imut = mut[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      ck = csix[k];
      ak = adisp[k];
      kmut = mut[k];

      constexpr bool incl = true;
      real xr = xi - xk;
      real yr = yi - yk;
      real zr = zi - zk;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real rr1 = REAL_RECIP(r);
         real e, de;
         real vlambda = 1;
         vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
         pair_disp<do_g, DTYP, 0, 1>(r, r2, rr1, scalea, aewald, ci, ai, ck, ak, vlambda, cut, off, e, de);
         if CONSTEXPR (do_e) {
            edtl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalea != 0 and e != 0)
                  ndtl += 1;
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de *= rr1;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            gxi += dedx;
            gyi += dedy;
            gzi += dedz;
            gxk -= dedx;
            gyk -= dedy;
            gzk -= dedz;
            if CONSTEXPR (do_v) {
               vdtlxx += floatTo<vbuf_prec>(xr * dedx);
               vdtlyx += floatTo<vbuf_prec>(yr * dedx);
               vdtlzx += floatTo<vbuf_prec>(zr * dedx);
               vdtlyy += floatTo<vbuf_prec>(yr * dedy);
               vdtlzy += floatTo<vbuf_prec>(zr * dedy);
               vdtlzz += floatTo<vbuf_prec>(zr * dedz);
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
      ci = csix[i];
      ai = adisp[i];
      imut = mut[i];
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;
      ck = csix[k];
      ak = adisp[k];
      kmut = mut[k];

      unsigned int dinfo0 = dinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dinfo0 & srcmask) == 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real rr1 = REAL_RECIP(r);
            real e, de;
            real vlambda = 1;
            vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
            pair_disp<do_g, DTYP, 1, 1>(r, r2, rr1, 1, aewald, ci, ai, ck, ak, vlambda, cut, off, e, de);
            if CONSTEXPR (do_e) {
               edtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     ndtl += 1;
               }
            }
            if CONSTEXPR (do_g) {
               real dedx, dedy, dedz;
               de *= rr1;
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               gxi += dedx;
               gyi += dedy;
               gzi += dedz;
               gxk -= dedx;
               gyk -= dedy;
               gzk -= dedz;
               if CONSTEXPR (do_v) {
                  vdtlxx += floatTo<vbuf_prec>(xr * dedx);
                  vdtlyx += floatTo<vbuf_prec>(yr * dedx);
                  vdtlzx += floatTo<vbuf_prec>(zr * dedx);
                  vdtlyy += floatTo<vbuf_prec>(yr * dedy);
                  vdtlzy += floatTo<vbuf_prec>(zr * dedy);
                  vdtlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ci = __shfl_sync(ALL_LANES, ci, ilane + 1);
         ai = __shfl_sync(ALL_LANES, ai, ilane + 1);
         imut = __shfl_sync(ALL_LANES, imut, ilane + 1);
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
      ci = csix[i];
      ai = adisp[i];
      imut = mut[i];
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;
      ck = csix[k];
      ak = adisp[k];
      kmut = mut[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         bool incl = atomk > 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real rr1 = REAL_RECIP(r);
            real e, de;
            real vlambda = 1;
            vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
            pair_disp<do_g, DTYP, 1, 1>(r, r2, rr1, 1, aewald, ci, ai, ck, ak, vlambda, cut, off, e, de);
            if CONSTEXPR (do_e) {
               edtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     ndtl += 1;
               }
            }
            if CONSTEXPR (do_g) {
               real dedx, dedy, dedz;
               de *= rr1;
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               gxi += dedx;
               gyi += dedy;
               gzi += dedz;
               gxk -= dedx;
               gyk -= dedy;
               gzk -= dedz;
               if CONSTEXPR (do_v) {
                  vdtlxx += floatTo<vbuf_prec>(xr * dedx);
                  vdtlyx += floatTo<vbuf_prec>(yr * dedx);
                  vdtlzx += floatTo<vbuf_prec>(zr * dedx);
                  vdtlyy += floatTo<vbuf_prec>(yr * dedy);
                  vdtlzy += floatTo<vbuf_prec>(zr * dedy);
                  vdtlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)

         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ci = __shfl_sync(ALL_LANES, ci, ilane + 1);
         ai = __shfl_sync(ALL_LANES, ai, ilane + 1);
         imut = __shfl_sync(ALL_LANES, imut, ilane + 1);
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
      atomic_add(ndtl, nd, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(edtl, ed, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vdtlxx, vdtlyx, vdtlzx, vdtlyy, vdtlzy, vdtlzz, vd, ithread);
   }
}
