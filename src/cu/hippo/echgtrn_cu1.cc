// ck.py Version 3.0.2
template <class Ver, Chgtrn CT>
__global__
void echgtrn_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nc, EnergyBuffer restrict ec,
   VirialBuffer restrict vc, grad_prec* restrict gx, grad_prec* restrict gy, grad_prec* restrict gz, real cut, real off,
   const unsigned* restrict minfo, int nexclude, const int (*restrict exclude)[2],
   const real (*restrict exclude_scale)[3], const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real* restrict chgct, real* restrict dmpct, real f)
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
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec ectl;
   if CONSTEXPR (do_e) {
      ectl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz;
   if CONSTEXPR (do_v) {
      vctlxx = 0;
      vctlyx = 0;
      vctlzx = 0;
      vctlyy = 0;
      vctlzy = 0;
      vctlzz = 0;
   }
   real xi, yi, zi, chgi, alphai;
   real xk, yk, zk, chgk, alphak;
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
      real scalea = exclude_scale[ii][0];

      xi = x[i];
      yi = y[i];
      zi = z[i];
      chgi = chgct[i];
      alphai = dmpct[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
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
         if CONSTEXPR (CT == Chgtrn::SEPARATE)
            pair_chgtrn<do_g>(r, cut, off, scalea, f, alphai, chgi, alphak, chgk, e, de);
         else if CONSTEXPR (CT == Chgtrn::COMBINED)
            pair_chgtrn_aplus<do_g>(r, cut, off, scalea, f, alphai, chgi, alphak, chgk, e, de);
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
      chgi = chgct[i];
      alphai = dmpct[i];
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;
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
            if CONSTEXPR (CT == Chgtrn::SEPARATE)
               pair_chgtrn<do_g>(r, cut, off, 1, f, alphai, chgi, alphak, chgk, e, de);
            else if CONSTEXPR (CT == Chgtrn::COMBINED)
               pair_chgtrn_aplus<do_g>(r, cut, off, 1, f, alphai, chgi, alphak, chgk, e, de);
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
      chgi = chgct[i];
      alphai = dmpct[i];
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;
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
            if CONSTEXPR (CT == Chgtrn::SEPARATE)
               pair_chgtrn<do_g>(r, cut, off, 1, f, alphai, chgi, alphak, chgk, e, de);
            else if CONSTEXPR (CT == Chgtrn::COMBINED)
               pair_chgtrn_aplus<do_g>(r, cut, off, 1, f, alphai, chgi, alphak, chgk, e, de);
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
