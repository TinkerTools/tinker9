#include "add.h"
#include "echarge.h"
#include "gpu_card.h"
#include "mdegv.h"
#include "mdpq.h"
#include "named_struct.h"
#include "nblist.h"
#include "pmestuf.h"
#include "seq_image.h"
#include "seq_pair_charge.h"
#include "seq_switch.h"
#include "switch.h"


TINKER_NAMESPACE_BEGIN
#define DEVICE_PTRS x, y, z, gx, gy, gz, pchg, nec, ec, vir_ec
template <class Ver, class ETYP>
void echarge_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   real f = electric / dielec;
   // real cut;
   real off, aewald;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      off = switch_off(switch_ewald);
      // cut = off;
      const auto& st = *epme_unit;
      aewald = st.aewald;
   } else if CONSTEXPR (eq<ETYP, NON_EWALD_TAPER>()) {
      off = switch_off(switch_charge);
      // cut = switch_cut(switch_charge);
   }
   const real off2 = off * off;
   const int maxnlist = clist_unit->maxnlst;
   const auto* nlst = clist_unit->nlst;
   const auto* lst = clist_unit->lst;
   auto bufsize = buffer_size();


   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(DEVICE_PTRS,nlst,lst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = pchg[i];
      MAYBE_UNUSED real gxi = 0, gyi = 0, gzi = 0;


      int nlsti = nlst[i];
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nlsti; ++kk) {
         int offset = (kk + i * n) & (bufsize - 1);
         int k = lst[i * maxnlist + kk];
         real xr = xi - x[k];
         real yr = yi - y[k];
         real zr = zi - z[k];
         real ck = pchg[k];


         real r2 = image2(xr, yr, zr);
         if (r2 <= off2) {
            real r = REAL_SQRT(r2);


            int ctl;
            real e, frcx, frcy, frcz, vxx, vxy, vxz, vyy, vyz, vzz;
            if CONSTEXPR (do_a) {
               ctl = 0;
            }
            if CONSTEXPR (do_e) {
               e = 0;
            }
            if CONSTEXPR (do_g) {
               frcx = 0;
               frcy = 0;
               frcz = 0;
            }
            if CONSTEXPR (do_v) {
               vxx = 0;
               vxy = 0;
               vxz = 0;
               vyy = 0;
               vyz = 0;
               vzz = 0;
            }


            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_charge<Ver, EWALD>(r, xr, yr, zr, 1, ci, ck, ebuffer, f,
                                       aewald, //
                                       frcx, frcy, frcz, ctl, e, vxx, vxy, vxz,
                                       vyy, vyz, vzz);
            }


            if CONSTEXPR (do_a)
               atomic_add(ctl, nec, offset);
            if CONSTEXPR (do_e)
               atomic_add(e, ec, offset);
            if CONSTEXPR (do_g) {
               gxi += frcx;
               gyi += frcy;
               gzi += frcz;
               atomic_add(-frcx, gx, k);
               atomic_add(-frcy, gy, k);
               atomic_add(-frcz, gz, k);
            }
            if CONSTEXPR (do_v)
               atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ec, offset);
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
      }
   } // end for (int i)


   #pragma acc parallel async deviceptr(DEVICE_PTRS, cexclude, cexclude_scale)
   #pragma acc loop independent
   for (int ii = 0; ii < ncexclude; ++ii) {
      int offset = ii & (bufsize - 1);


      int i = cexclude[ii][0];
      int k = cexclude[ii][1];
      real cscale = cexclude_scale[ii];


      real ci = pchg[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real ck = pchg[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         int ctl;
         real e, frcx, frcy, frcz, vxx, vxy, vxz, vyy, vyz, vzz;
         if CONSTEXPR (do_a) {
            ctl = 0;
         }
         if CONSTEXPR (do_e) {
            e = 0;
         }
         if CONSTEXPR (do_g) {
            frcx = 0;
            frcy = 0;
            frcz = 0;
         }
         if CONSTEXPR (do_v) {
            vxx = 0;
            vxy = 0;
            vxz = 0;
            vyy = 0;
            vyz = 0;
            vzz = 0;
         }


         real r = REAL_SQRT(r2);
         if CONSTEXPR (eq<ETYP, EWALD>()) {
            pair_charge<Ver, NON_EWALD>(
               r, xr, yr, zr, cscale, ci, ck, ebuffer, f, 0, //
               frcx, frcy, frcz, ctl, e, vxx, vxy, vxz, vyy, vyz, vzz);
         }


         if (e != 0) {
            if CONSTEXPR (do_a)
               atomic_add(ctl, nec, offset);
            if CONSTEXPR (do_e)
               atomic_add(e, ec, offset);
         }
         if CONSTEXPR (do_g) {
            atomic_add(frcx, gx, i);
            atomic_add(frcy, gy, i);
            atomic_add(frcz, gz, i);
            atomic_add(-frcx, gx, k);
            atomic_add(-frcy, gy, k);
            atomic_add(-frcz, gz, k);
         }
         if CONSTEXPR (do_v)
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ec, offset);
      } // end if (include)
   }
}


void echarge_ewald_real_acc(int vers)
{
   if (vers == calc::v0)
      echarge_acc1<calc::V0, EWALD>();
   else if (vers == calc::v1)
      echarge_acc1<calc::V1, EWALD>();
   else if (vers == calc::v3)
      echarge_acc1<calc::V3, EWALD>();
   else if (vers == calc::v4)
      echarge_acc1<calc::V4, EWALD>();
   else if (vers == calc::v5)
      echarge_acc1<calc::V5, EWALD>();
   else if (vers == calc::v6)
      echarge_acc1<calc::V6, EWALD>();
}
TINKER_NAMESPACE_END
