
template <class Ver>
__global__
void echgtrn_cu1(
   int n, TINKER_IMAGE_PARAMS, count_buffer restrict nc,
   energy_buffer restrict ec, virial_buffer restrict vc, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, real cut, real off,
   const unsigned* restrict minfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
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
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec ectl;
   if CONSTEXPR (do_e) {
      ectl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz;
   if CONSTEXPR (do_v) {
      vctlxx = 0;
      vctlyx = 0;
      vctlzx = 0;
      vctlyy = 0;
      vctlzy = 0;
      vctlzz = 0;
   }


   real shxi;
   real shyi;
   real shzi;
   real xk;
   real yk;
   real zk;
   real shgxi;
   real shgyi;
   real shgzi;
   real gxk;
   real gyk;
   real gzk;
   real shchgi;
   real shalphai;
   real chgk;
   real alphak;


   //* /
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         shgxi = 0;
         shgyi = 0;
         shgzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
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
      real chgi = chgct[shi];
      real alphai = dmpct[shi];
      chgk = chgct[k];
      alphak = dmpct[k];


      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;

      MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);

         MAYBE_UNUSED e_prec e, de;
         pair_chgtrn<do_g>(r, scalea, f, alphai, chgi, alphak, chgk, e, de);


         if (r2 > cut * cut) {
            real taper, dtaper;
            switch_taper5<do_g>(r, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e *= taper;
         }


         if CONSTEXPR (do_a)
            if (e != 0)
               nctl += 1;
         if CONSTEXPR (do_e)
            ectl += cvt_to<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(r);
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;

            shgxi += dedx;
            shgyi += dedy;
            shgzi += dedz;
            gxk -= dedx;
            gyk -= dedy;
            gzk -= dedz;

            if CONSTEXPR (do_v) {
               vctlxx += cvt_to<vbuf_prec>(xr * dedx);
               vctlyx += cvt_to<vbuf_prec>(yr * dedx);
               vctlzx += cvt_to<vbuf_prec>(zr * dedx);
               vctlyy += cvt_to<vbuf_prec>(yr * dedy);
               vctlzy += cvt_to<vbuf_prec>(zr * dedy);
               vctlzz += cvt_to<vbuf_prec>(zr * dedz);
            }
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shgxi, gx, shi);
         atomic_add(shgyi, gy, shi);
         atomic_add(shgzi, gz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi = 0;
         shgyi = 0;
         shgzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
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
      shxi = sorted[shatomi].x;
      shyi = sorted[shatomi].y;
      shzi = sorted[shatomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;


      shchgi = chgct[shi];
      shalphai = dmpct[shi];
      chgk = chgct[k];
      alphak = dmpct[k];


      unsigned int minfo0 = minfo[iw * WARP_SIZE + ilane];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int srcmask = 1 << srclane;

         int iid = shiid;
         real xi = shxi;
         real yi = shyi;
         real zi = shzi;
         real chgi = shchgi;
         real alphai = shalphai;


         bool incl = iid < kid and kid < n;
         incl = incl and (minfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);

            MAYBE_UNUSED e_prec e, de;
            pair_chgtrn<do_g>(r, scalea, f, alphai, chgi, alphak, chgk, e, de);


            if (r2 > cut * cut) {
               real taper, dtaper;
               switch_taper5<do_g>(r, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e *= taper;
            }


            if CONSTEXPR (do_a)
               if (e != 0)
                  nctl += 1;
            if CONSTEXPR (do_e)
               ectl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(r);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;

               shgxi += dedx;
               shgyi += dedy;
               shgzi += dedz;
               gxk -= dedx;
               gyk -= dedy;
               gzk -= dedz;

               if CONSTEXPR (do_v) {
                  vctlxx += cvt_to<vbuf_prec>(xr * dedx);
                  vctlyx += cvt_to<vbuf_prec>(yr * dedx);
                  vctlzx += cvt_to<vbuf_prec>(zr * dedx);
                  vctlyy += cvt_to<vbuf_prec>(yr * dedy);
                  vctlzy += cvt_to<vbuf_prec>(zr * dedy);
                  vctlzz += cvt_to<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
         shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
         shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
         shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
         shchgi = __shfl_sync(ALL_LANES, shchgi, ilane + 1);
         shalphai = __shfl_sync(ALL_LANES, shalphai, ilane + 1);
         if CONSTEXPR (do_g) {
            shgxi = __shfl_sync(ALL_LANES, shgxi, ilane + 1);
            shgyi = __shfl_sync(ALL_LANES, shgyi, ilane + 1);
            shgzi = __shfl_sync(ALL_LANES, shgzi, ilane + 1);
         }
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi, gx, shi);
         atomic_add(shgyi, gy, shi);
         atomic_add(shgzi, gz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi = 0;
         shgyi = 0;
         shgzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }


      int ty = iak[iw];
      int shatomi = ty * WARP_SIZE + ilane;
      int shi = sorted[shatomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      shxi = sorted[shatomi].x;
      shyi = sorted[shatomi].y;
      shzi = sorted[shatomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;


      shchgi = chgct[shi];
      shalphai = dmpct[shi];
      chgk = chgct[k];
      alphak = dmpct[k];


      for (int j = 0; j < WARP_SIZE; ++j) {

         real xi = shxi;
         real yi = shyi;
         real zi = shzi;
         real chgi = shchgi;
         real alphai = shalphai;


         bool incl = atomk > 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);

            MAYBE_UNUSED e_prec e, de;
            pair_chgtrn<do_g>(r, scalea, f, alphai, chgi, alphak, chgk, e, de);


            if (r2 > cut * cut) {
               real taper, dtaper;
               switch_taper5<do_g>(r, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e *= taper;
            }


            if CONSTEXPR (do_a)
               if (e != 0)
                  nctl += 1;
            if CONSTEXPR (do_e)
               ectl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(r);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;

               shgxi += dedx;
               shgyi += dedy;
               shgzi += dedz;
               gxk -= dedx;
               gyk -= dedy;
               gzk -= dedz;

               if CONSTEXPR (do_v) {
                  vctlxx += cvt_to<vbuf_prec>(xr * dedx);
                  vctlyx += cvt_to<vbuf_prec>(yr * dedx);
                  vctlzx += cvt_to<vbuf_prec>(zr * dedx);
                  vctlyy += cvt_to<vbuf_prec>(yr * dedy);
                  vctlzy += cvt_to<vbuf_prec>(zr * dedy);
                  vctlzz += cvt_to<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)


         shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
         shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
         shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
         shchgi = __shfl_sync(ALL_LANES, shchgi, ilane + 1);
         shalphai = __shfl_sync(ALL_LANES, shalphai, ilane + 1);
         if CONSTEXPR (do_g) {
            shgxi = __shfl_sync(ALL_LANES, shgxi, ilane + 1);
            shgyi = __shfl_sync(ALL_LANES, shgyi, ilane + 1);
            shgzi = __shfl_sync(ALL_LANES, shgzi, ilane + 1);
         }
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi, gx, shi);
         atomic_add(shgyi, gy, shi);
         atomic_add(shgzi, gz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }
   // */


   if CONSTEXPR (do_a) {
      atomic_add(nctl, nc, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(ectl, ec, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz, vc, ithread);
   }
} // generated by ComplexKernelBuilder (ck.py) 1.5.1
