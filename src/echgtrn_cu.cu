#include "add.h"
#include "echgtrn.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "mod.chgpot.h"
#include "seq_pair_chgtrn.h"
#include "seq_switch.h"
#include "switch.h"
#include <cassert>


namespace tinker {
#define HIPPO_CHGTRN_PARA                                                      \
   size_t bufsize, count_buffer restrict nct, energy_buffer restrict ect,      \
      virial_buffer restrict vir_ect, grad_prec *restrict gx, grad_prec *gy,   \
      grad_prec *gz, TINKER_IMAGE_PARAMS, real *restrict chgct,                \
      real *restrict dmpct, real f, real cut, real off


template <class Ver>
__global__
void echgtrn_cu1(HIPPO_CHGTRN_PARA, int n,
                 const Spatial::SortedAtom* restrict sorted, int niak,
                 const int* restrict iak, const int* restrict lst)
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
      real chgi = chgct[i];
      real alphai = dmpct[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      real shchgk = chgct[shk];
      real shalphak = dmpct[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = __shfl_sync(ALL_LANES, shx, srclane) - xi;
         real yr = __shfl_sync(ALL_LANES, shy, srclane) - yi;
         real zr = __shfl_sync(ALL_LANES, shz, srclane) - zi;
         real alphak = __shfl_sync(ALL_LANES, shalphak, srclane);
         real chgk = __shfl_sync(ALL_LANES, shchgk, srclane);


         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         real r2 = image2(xr, yr, zr);
         if (atomi < atomk && r2 <= off2) {
            real r = REAL_SQRT(r2);


            MAYBE_UNUSED e_prec e, de;
            pair_chgtrn<do_g>(r, 1, f, alphai, chgi, alphak, chgk, e, de);


            if (r2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(r, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e *= taper;
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
            gxi -= dedx;
            gyi -= dedy;
            gzi -= dedz;
            gxk += __shfl_sync(ALL_LANES, dedx, dstlane);
            gyk += __shfl_sync(ALL_LANES, dedy, dstlane);
            gzk += __shfl_sync(ALL_LANES, dedz, dstlane);
         }
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, nct, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, ect, offset);
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, shk);
         atomic_add(gyk, gy, shk);
         atomic_add(gzk, gz, shk);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_ect, offset);
   } // end for (iw)
}


template <class Ver>
__global__
void echgtrn_cu2(HIPPO_CHGTRN_PARA, const real* x, const real* y, const real* z,
                 int nmexclude, int (*restrict mexclude)[2],
                 real* restrict mexclude_scale)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nmexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = mexclude[ii][0];
      int k = mexclude[ii][1];
      real mscale = mexclude_scale[ii];


      real alphai = dmpct[i];
      real chgi = chgct[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real alphak = dmpct[k];
      real chgk = chgct[k];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real r = REAL_SQRT(r2);


         MAYBE_UNUSED e_prec e, de;
         pair_chgtrn<do_g>(r, mscale, f, alphai, chgi, alphak, chgk, e, de);
         if (r2 > cut2) {
            real taper, dtaper;
            switch_taper5<do_g>(r, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e = e * taper;
         }


         if CONSTEXPR (do_a)
            if (mscale == -1 && e != 0)
               atomic_add(-1, nct, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, ect, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(r);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(-dedx, gx, i);
            atomic_add(-dedy, gy, i);
            atomic_add(-dedz, gz, i);
            atomic_add(dedx, gx, k);
            atomic_add(dedy, gy, k);
            atomic_add(dedz, gz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ect, offset);
            }
         }
      }
   }
}


template <class Ver>
void echgtrn_cu3()
{
   const auto& st = *mspatial_unit;
   real cut = switch_cut(switch_chgtrn);
   real off = switch_off(switch_chgtrn);
   real f = electric / dielec;


   auto bufsize = buffer_size();


   assert(ctrntyp == chgtrn_t::SEPARATE);


   auto ker1 = echgtrn_cu1<Ver>;
   if (st.niak > 0)
      launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                 bufsize, nct, ect, vir_ect, dectx, decty, dectz,
                 TINKER_IMAGE_ARGS, chgct, dmpct, f, cut, off, //
                 n, st.sorted, st.niak, st.iak, st.lst);


   auto ker2 = echgtrn_cu2<Ver>;
   if (nmexclude > 0)
      launch_k1s(nonblk, nmexclude, ker2, //
                 bufsize, nct, ect, vir_ect, dectx, decty, dectz,
                 TINKER_IMAGE_ARGS, chgct, dmpct, f, cut, off, x, y, z,
                 nmexclude, mexclude, mexclude_scale);
}


void echgtrn_cu(int vers)
{
   if (vers == calc::v0)
      echgtrn_cu3<calc::V0>();
   else if (vers == calc::v1)
      echgtrn_cu3<calc::V1>();
   else if (vers == calc::v3)
      echgtrn_cu3<calc::V3>();
   else if (vers == calc::v4)
      echgtrn_cu3<calc::V4>();
   else if (vers == calc::v5)
      echgtrn_cu3<calc::V5>();
   else if (vers == calc::v6)
      echgtrn_cu3<calc::V6>();
}
}
