#include "add.h"
#include "ehippo.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "mdpq.h"
#include "named_struct.h"
#include "seq_switch.h"
#include "switch.h"


namespace tinker {
template <class Ver>
__device__
void pair_hippo(real r2, real3 dr, real mscale, real alphai, real chgi,
                real alphak, real chgk, real cut, real off, real f,
                real3& restrict frci, real3& restrict frck, real& restrict etl,
                real& restrict vtlxx, real& restrict vtlxy,
                real& restrict vtlxz, real& restrict vtlyy,
                real& restrict vtlyz, real& restrict vtlzz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   // HIPPO charge transfer uses 'separate' type only.
   real r = REAL_SQRT(r2);
   real expi = REAL_EXP(-alphai * r);
   real expk = REAL_EXP(-alphak * r);


   real e, de;
   e = -chgi * expk - chgk * expi;
   if CONSTEXPR (do_g) {
      de = chgi * expk * alphak + chgk * expi * alphai;
   }


   if (r > cut) {
      real taper, dtaper;
      switch_taper5<do_g>(r, cut, off, taper, dtaper);
      if CONSTEXPR (do_g)
         de = e * dtaper + de * taper;
      if CONSTEXPR (do_e)
         e = e * taper;
   }


   if CONSTEXPR (do_e) {
      e *= f * mscale;
      etl = etl + e;
   }
   real frcx, frcy, frcz;
   if CONSTEXPR (do_g) {
      de *= f * mscale;
      real rr1 = REAL_RECIP(r);
      frcx = de * dr.x * rr1;
      frcy = de * dr.y * rr1;
      frcz = de * dr.z * rr1;
      frci.x = frci.x - frcx;
      frci.y = frci.y - frcy;
      frci.z = frci.z - frcz;
      frck.x = frck.x + frcx;
      frck.y = frck.y + frcy;
      frck.z = frck.z + frcz;
   }
   if CONSTEXPR (do_v) {
      vtlxx += dr.x * frcx;
      vtlxy += dr.y * frcx;
      vtlxz += dr.z * frcx;
      vtlyy += dr.y * frcy;
      vtlyz += dr.z * frcy;
      vtlzz += dr.z * frcz;
   }
}


//====================================================================//


#define HIPPO_PARA                                                             \
   size_t bufsize, count_buffer restrict nct, energy_buffer restrict ect,      \
      virial_buffer restrict vir_ect, grad_prec *restrict gx, grad_prec *gy,   \
      grad_prec *gz, TINKER_IMAGE_PARAMS, real *restrict chgct,                \
      real *restrict dmpct, real cut, real off


template <class Ver>
__global__
void ehippo_cu1(HIPPO_PARA, int n, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst)
{}


template <class Ver>
__global__
void ehippo_cu2(HIPPO_PARA, const real* x, const real* y, const real* z,
                int nmexclude, int (*restrict mexlcude)[2],
                real* restrict mexclude_scale)
{}


template <class Ver>
void ehippo_cu3()
{
   const auto& st = *mspatial_unit;
   real cut = switch_cut(switch_chgtrn);
   real off = switch_off(switch_chgtrn);


   auto bufsize = buffer_size();


   auto ker1 = ehippo_cu1<Ver>;
   launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
              bufsize, nct, ect, vir_ect, dectx, decty, dectz,
              TINKER_IMAGE_ARGS, chgct, dmpct, cut, off, //
              n, st.sorted, st.niak, st.iak, st.lst);


   auto ker2 = ehippo_cu2<Ver>;
   launch_k1s(nonblk, nmexclude, ker2, //
              bufsize, nct, ect, vir_ect, dectx, decty, dectz,
              TINKER_IMAGE_ARGS, chgct, dmpct, cut, off, x, y, z, nmexclude,
              mexclude, mexclude_scale);
}


void ehippo_cu(int vers)
{
   if (vers == calc::v0)
      ehippo_cu3<calc::V0>();
   else if (vers == calc::v1)
      ehippo_cu3<calc::V1>();
   // else if (vers == calc::v3)
   //    ehippo_cu3<calc::V3>();
   else if (vers == calc::v4)
      ehippo_cu3<calc::V4>();
   else if (vers == calc::v5)
      ehippo_cu3<calc::V5>();
   else if (vers == calc::v6)
      ehippo_cu3<calc::V6>();
}
}
