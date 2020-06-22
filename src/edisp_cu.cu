#include "add.h"
#include "edisp.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pmestuf.h"
#include "seq_switch.h"
#include "switch.h"


namespace tinker {
template <bool DO_G, class DTYP>
__device__
void pair_disp(real r, real r2, real dspscale, real aewald, real ci, real ai,
               real ck, real ak, e_prec& restrict e, e_prec& restrict de)
{
   real rr1 = REAL_RECIP(r);
   real rr2 = REAL_RECIP(r2);
   real rr6 = rr2 * rr2 * rr2;


   real di = ai * r;
   real dk = ak * r;
   real expi = REAL_EXP(-di);
   real expk = REAL_EXP(-dk);
   real di2 = di * di;
   real di3 = di * di2;
   real damp3 = 1, damp5 = 1, ddamp = 0;
   if (ai != ak) {
      real ai2 = ai * ai;
      real ak2 = ak * ak;
      real dk2 = dk * dk;
      real dk3 = dk * dk2;
      real ti = ak2 * REAL_RECIP(ak2 - ai2);
      real tk = ai2 * REAL_RECIP(ai2 - ak2);
      real ti2 = ti * ti;
      real tk2 = tk * tk;
      damp3 = 1 - ti2 * (1 + di + 0.5f * di2) * expi -
         tk2 * (1 + dk + 0.5f * dk2) * expk - 2 * ti2 * tk * (1 + di) * expi -
         2 * tk2 * ti * (1 + dk) * expk;
      damp5 = 1 - ti2 * (1 + di + 0.5f * di2 + di3 / 6) * expi -
         tk2 * (1 + dk + 0.5f * dk2 + dk3 / 6) * expk -
         2 * ti2 * tk * (1 + di + di2 / 3) * expi -
         2 * tk2 * ti * (1 + dk + dk2 / 3) * expk;
      if CONSTEXPR (DO_G)
         ddamp = 0.25f * di2 * ti2 * ai * expi * (r * ai + 4 * tk - 1) +
            0.25f * dk2 * tk2 * ak * expk * (r * ak + 4 * ti - 1);
   }
   if (ai == ak) {
      real di4 = di2 * di2;
      real di5 = di2 * di3;
      damp3 = 1 - (1 + di + 0.5f * di2 + 7 * di3 / 48 + di4 / 48) * expi;
      damp5 = 1 - (1 + di + 0.5f * di2 + di3 / 6 + di4 / 24 + di5 / 144) * expi;
      if CONSTEXPR (DO_G)
         ddamp = ai * expi * (di5 - 3 * di3 - 3 * di2) / 96;
   }
   real damp = 1.5f * damp5 - 0.5f * damp3;


   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      real ralpha2 = r2 * aewald * aewald;
      real term = 1 + ralpha2 + 0.5f * ralpha2 * ralpha2;
      real expterm = REAL_EXP(-ralpha2);
      real expa = expterm * term;
      e = -ci * ck * rr6 * (expa + damp * damp - 1);
      if CONSTEXPR (DO_G) {
         real rterm = -ralpha2 * ralpha2 * ralpha2 * rr1 * expterm;
         de = -6 * e * rr1 - ci * ck * rr6 * (rterm + 2 * damp * ddamp);
      }
   } else if CONSTEXPR (eq<DTYP, NON_EWALD>() || eq<DTYP, NON_EWALD_TAPER>()) {
      e = -ci * ck * rr6;
      if CONSTEXPR (DO_G) {
         de = -6 * e * rr1;
         de = dspscale * (de * damp * damp + 2 * e * damp * ddamp);
      }
      e = dspscale * (e * damp * damp);
   }
}


//====================================================================//


#define EDISP_ARGS                                                             \
   size_t bufsize, count_buffer restrict ndisp, energy_buffer restrict edsp,   \
      virial_buffer restrict vir_edsp, grad_prec *restrict gx,                 \
      grad_prec *restrict gy, grad_prec *restrict gz, TINKER_IMAGE_PARAMS,     \
      real cut, real off, const real *restrict csix,                           \
      const real *restrict adisp


template <class Ver, class DTYP>
__global__
void edisp_cu1(EDISP_ARGS, const Spatial::SortedAtom* restrict sorted, int niak,
               const int* restrict iak, const int* restrict lst, int n,
               real aewald)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   MAYBE_UNUSED int ctl;
   MAYBE_UNUSED e_prec etl;
   MAYBE_UNUSED v_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
   MAYBE_UNUSED real gxi, gyi, gzi, gxk, gyk, gzk;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_a)
         ctl = 0;
      if CONSTEXPR (do_e)
         etl = 0;
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlyx = 0;
         vtlzx = 0;
         vtlyy = 0;
         vtlzy = 0;
         vtlzz = 0;
      }


      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      real ci = csix[i];
      real ai = adisp[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      real shck = csix[shk];
      real shak = adisp[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = xi - __shfl_sync(ALL_LANES, shx, srclane);
         real yr = yi - __shfl_sync(ALL_LANES, shy, srclane);
         real zr = zi - __shfl_sync(ALL_LANES, shz, srclane);
         real ck = __shfl_sync(ALL_LANES, shck, srclane);
         real ak = __shfl_sync(ALL_LANES, shak, srclane);


         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         real r2 = image2(xr, yr, zr);
         if (atomi < atomk && r2 <= off2) {
            real r = REAL_SQRT(r2);


            MAYBE_UNUSED e_prec e, de;
            if CONSTEXPR (eq<DTYP, DEWALD>()) {
               pair_disp<do_g, DEWALD>(r, r2, 1, aewald, ci, ai, ck, ak, e, de);
            } else if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
               pair_disp<do_g, NON_EWALD_TAPER>(r, r2, 1, 0, ci, ai, ck, ak, e,
                                                de);
               if (r2 > cut2) {
                  real taper, dtaper;
                  switch_taper5<do_g>(r, cut, off, taper, dtaper);
                  if CONSTEXPR (do_g)
                     de = e * dtaper + de * taper;
                  if CONSTEXPR (do_e)
                     e *= taper;
               }
            }


            if CONSTEXPR (do_a)
               if (e != 0)
                  ctl += 1;
            if CONSTEXPR (do_e)
               etl += e;
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(r);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               if CONSTEXPR (do_v) {
                  vtlxx += xr * dedx;
                  vtlyx += yr * dedx;
                  vtlzx += zr * dedx;
                  vtlyy += yr * dedy;
                  vtlzy += zr * dedy;
                  vtlzz += zr * dedz;
               }
            }
         } // end if (include)


         if CONSTEXPR (do_g) {
            int dstlane = (ilane + WARP_SIZE - j) & (WARP_SIZE - 1);
            gxi += dedx;
            gyi += dedy;
            gzi += dedz;
            gxk -= __shfl_sync(ALL_LANES, dedx, dstlane);
            gyk -= __shfl_sync(ALL_LANES, dedy, dstlane);
            gzk -= __shfl_sync(ALL_LANES, dedz, dstlane);
         }
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, ndisp, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, edsp, offset);
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, shk);
         atomic_add(gyk, gy, shk);
         atomic_add(gzk, gz, shk);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_edsp, offset);
   } // end for (iw)
}


template <class Ver, class DTYP>
__global__
void edisp_cu2(EDISP_ARGS, const real* restrict x, const real* restrict y,
               const real* restrict z, int ndspexclude,
               const int (*restrict dspexclude)[2],
               const real* restrict dspexclude_scale)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndspexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = dspexclude[ii][0];
      int k = dspexclude[ii][1];
      real dspscale = dspexclude_scale[ii];


      real ci = csix[i];
      real ai = adisp[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real ck = csix[k];
      real ak = adisp[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real r = REAL_SQRT(r2);


         MAYBE_UNUSED e_prec e, de;
         pair_disp<do_g, NON_EWALD>(r, r2, dspscale, 0, ci, ai, ck, ak, e, de);
         if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
            if (r2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(r, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e = e * taper;
            }
         }


         if CONSTEXPR (do_a)
            if (dspscale == -1 && e != 0)
               atomic_add(-1, ndisp, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, edsp, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(r);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, gx, i);
            atomic_add(dedy, gy, i);
            atomic_add(dedz, gz, i);
            atomic_add(-dedx, gx, k);
            atomic_add(-dedy, gy, k);
            atomic_add(-dedz, gz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_edsp, offset);
            }
         }
      }
   }
}


template <class Ver, class DTYP>
void edisp_cu()
{
   const auto& st = *dspspatial_unit;
   real cut, off;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      off = switch_off(switch_dewald);
      // cut = off; // not used
   } else {
      off = switch_off(switch_disp);
      cut = switch_cut(switch_disp);
   }
   size_t bufsize = buffer_size();


   real aewald = 0;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      PMEUnit pu = dpme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      if (st.niak > 0) {
         auto ker1 = edisp_cu1<Ver, DEWALD>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, 0, off, csix, adisp, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald);
      }
      if (ndspexclude > 0) {
         auto ker2 = edisp_cu2<Ver, NON_EWALD>;
         launch_k1s(nonblk, ndspexclude, ker2, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, 0, off, csix, adisp, //
                    x, y, z, ndspexclude, dspexclude, dspexclude_scale);
      }
   } else if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
      if (st.niak > 0) {
         auto ker1 = edisp_cu1<Ver, NON_EWALD_TAPER>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, cut, off, csix, adisp, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald);
      }
      if (ndspexclude > 0) {
         auto ker2 = edisp_cu2<Ver, NON_EWALD_TAPER>;
         launch_k1s(nonblk, ndspexclude, ker2, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, cut, off, csix, adisp, //
                    x, y, z, ndspexclude, dspexclude, dspexclude_scale);
      }
   }
}


void edisp_ewald_real_cu(int vers)
{
   if (vers == calc::v0)
      edisp_cu<calc::V0, DEWALD>();
   else if (vers == calc::v1)
      edisp_cu<calc::V1, DEWALD>();
   else if (vers == calc::v3)
      edisp_cu<calc::V3, DEWALD>();
   else if (vers == calc::v4)
      edisp_cu<calc::V4, DEWALD>();
   else if (vers == calc::v5)
      edisp_cu<calc::V5, DEWALD>();
   else if (vers == calc::v6)
      edisp_cu<calc::V6, DEWALD>();
}


void edisp_nonewald_cu(int vers)
{
   if (vers == calc::v0)
      edisp_cu<calc::V0, NON_EWALD_TAPER>();
   else if (vers == calc::v1)
      edisp_cu<calc::V1, NON_EWALD_TAPER>();
   else if (vers == calc::v3)
      edisp_cu<calc::V3, NON_EWALD_TAPER>();
   else if (vers == calc::v4)
      edisp_cu<calc::V4, NON_EWALD_TAPER>();
   else if (vers == calc::v5)
      edisp_cu<calc::V5, NON_EWALD_TAPER>();
   else if (vers == calc::v6)
      edisp_cu<calc::V6, NON_EWALD_TAPER>();
}
}
