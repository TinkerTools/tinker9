#include "ff/atom.h"
#include "ff/evdw.h"
#include "ff/image.h"
#include "ff/nblist.h"
#include "ff/switch.h"
#include "math/switch.h"
#include "seq/add.h"
#include "seq/pair_lj.h"
#include "tool/externfunc.h"
#include "tool/gpucard.h"

namespace tinker {
#define DEVICE_PTRS x, y, z, devx, devy, devz, jvdw, radmin, epsilon, mut, nev, ev, vir_ev
template <class Ver>
static void elj_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   const real cut = switchCut(Switch::VDW);
   const real off = switchOff(Switch::VDW);
   const real cut2 = cut * cut;
   const real off2 = off * off;
   const int maxnlist = clist_unit->maxnlst;
   const auto* nlst = clist_unit->nlst;
   const auto* lst = clist_unit->lst;

   auto bufsize = bufferSize();

   MAYBE_UNUSED int GRID_DIM = gpuGridSize(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(DEVICE_PTRS,nlst,lst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      int it = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real imut = mut[i];
      MAYBE_UNUSED grad_prec gxi = 0, gyi = 0, gzi = 0;

      int nvlsti = nlst[i];
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nvlsti; ++kk) {
         int offset = (kk + i * n) & (bufsize - 1);
         int k = lst[i * maxnlist + kk];
         int kt = jvdw[k];
         real xr = xi - x[k];
         real yr = yi - y[k];
         real zr = zi - z[k];
         int kmut = mut[k];

         real rik2 = image2(xr, yr, zr);
         if (rik2 <= off2) {
            real rik = REAL_SQRT(rik2);
            real rv = radmin[it * njvdw + kt];
            real eps = epsilon[it * njvdw + kt];

            MAYBE_UNUSED real e, de;
            real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
            pair_lj_v1<do_g, true>(rik, vlambda, rv, eps, 1, e, de);

            if (rik2 > cut2) {
               real taper, dtaper;
               switchTaper5<do_g>(rik, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e = e * taper;
            }

            if CONSTEXPR (do_a)
               // vscale is always 1, exclude e == 0
               if (e != 0)
                  atomic_add(1, nev, offset);
            if CONSTEXPR (do_e)
               atomic_add(e, ev, offset);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(rik);
               real dedx = de * xr;
               real dedy = de * yr;
               real dedz = de * zr;
               gxi += floatTo<grad_prec>(dedx);
               gyi += floatTo<grad_prec>(dedy);
               gzi += floatTo<grad_prec>(dedz);
               atomic_add(-dedx, devx, k);
               atomic_add(-dedy, devy, k);
               atomic_add(-dedz, devz, k);
               if CONSTEXPR (do_v) {
                  real vxx = xr * dedx;
                  real vyx = yr * dedx;
                  real vzx = zr * dedx;
                  real vyy = yr * dedy;
                  real vzy = zr * dedy;
                  real vzz = zr * dedz;
                  atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
               }
            } // end if (do_g)
         }
      } // end for (int kk)

      if CONSTEXPR (do_g) {
         atomic_add(gxi, devx, i);
         atomic_add(gyi, devy, i);
         atomic_add(gzi, devz, i);
      }
   } // enf for (int i)

   #pragma acc parallel async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(DEVICE_PTRS,vexclude,vexclude_scale)
   #pragma acc loop independent
   for (int ii = 0; ii < nvexclude; ++ii) {
      int offset = ii & (bufsize - 1);

      int i = vexclude[ii][0];
      int k = vexclude[ii][1];
      real vscale = vexclude_scale[ii] - 1;

      int it = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      int imut = mut[i];

      int kt = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];
      int kmut = mut[k];

      real rik2 = image2(xr, yr, zr);
      if (rik2 <= off2) {
         real rik = REAL_SQRT(rik2);
         real rv = radmin[it * njvdw + kt];
         real eps = epsilon[it * njvdw + kt];

         MAYBE_UNUSED real e, de;
         real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
         pair_lj_v1<do_g, true>(rik, vlambda, rv, eps, vscale, e, de);

         if (rik2 > cut2) {
            real taper, dtaper;
            switchTaper5<do_g>(rik, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e = e * taper;
         }

         if CONSTEXPR (do_a)
            if (vscale == -1 && e != 0)
               atomic_add(-1, nev, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, ev, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(rik);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, devx, i);
            atomic_add(dedy, devy, i);
            atomic_add(dedz, devz, i);
            atomic_add(-dedx, devx, k);
            atomic_add(-dedy, devy, k);
            atomic_add(-dedz, devz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
            }
         } // end if (do_g)
      }
   } // end for (int ii)

   #pragma acc parallel async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(DEVICE_PTRS,vdw14ik,radmin4,epsilon4)
   #pragma acc loop independent
   for (int ii = 0; ii < nvdw14; ++ii) {
      int offset = ii & (bufsize - 1);

      int i = vdw14ik[ii][0];
      int k = vdw14ik[ii][1];

      int it = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      int imut = mut[i];

      int kt = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];
      int kmut = mut[k];

      real rik2 = image2(xr, yr, zr);
      if (rik2 <= off2) {
         real rik = REAL_SQRT(rik2);
         real rv = radmin[it * njvdw + kt];
         real eps = epsilon[it * njvdw + kt];
         real rv4 = radmin4[it * njvdw + kt];
         real eps4 = epsilon4[it * njvdw + kt];

         MAYBE_UNUSED real e, de, e4, de4;
         real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
         pair_lj_v1<do_g, true>(rik, vlambda, rv, eps, v4scale, e, de);
         pair_lj_v1<do_g, true>(rik, vlambda, rv4, eps4, v4scale, e4, de4);
         e = e4 - e;
         if CONSTEXPR (do_g)
            de = de4 - de;

         if (rik2 > cut2) {
            real taper, dtaper;
            switchTaper5<do_g>(rik, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e = e * taper;
         }

         // if CONSTEXPR (do_a) {}
         if CONSTEXPR (do_e)
            atomic_add(e, ev, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(rik);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, devx, i);
            atomic_add(dedy, devy, i);
            atomic_add(dedz, devz, i);
            atomic_add(-dedx, devx, k);
            atomic_add(-dedy, devy, k);
            atomic_add(-dedz, devz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
            }
         } // end if (do_g)
      }
   } // end for (int ii)
}

void elj_acc(int vers)
{
   if (vers == calc::v0)
      elj_acc1<calc::V0>();
   else if (vers == calc::v1)
      elj_acc1<calc::V1>();
   else if (vers == calc::v3)
      elj_acc1<calc::V3>();
   else if (vers == calc::v4)
      elj_acc1<calc::V4>();
   else if (vers == calc::v5)
      elj_acc1<calc::V5>();
   else if (vers == calc::v6)
      elj_acc1<calc::V6>();
}

void ebuck_acc(int)
{
   TINKER_F1CALL_ERROR__(ebuck, acc);
}

void emm3hb_acc(int)
{
   TINKER_F1CALL_ERROR__(emm3hb, acc);
}

void egauss_acc(int)
{
   TINKER_F1CALL_ERROR__(egauss, acc);
}
}
