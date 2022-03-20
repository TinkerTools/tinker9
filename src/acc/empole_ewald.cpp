#include "add.h"
#include "amoeba/empole.h"
#include "glob/nblist.h"
#include "image.h"
#include "md.h"
#include "pmestuf.h"
#include "seq/pair_mpole.h"
#include "switch.h"
#include "tool/gpucard.h"

namespace tinker {
#define DEVICE_PTRS x, y, z, demx, demy, demz, rpole, nem, em, vir_em, trqx, trqy, trqz
template <class Ver>
void empole_ewald_real_self_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   const real f = electric / dielec;

   const real off = switch_off(switch_ewald);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   auto bufsize = buffer_size();
   PairMPoleGrad pgrad;

   const PMEUnit pu = epme_unit;
   const real aewald = pu->aewald;
   const real aewald_sq_2 = 2 * aewald * aewald;
   const real fterm = -f * aewald * 0.5f * (real)(M_2_SQRTPI);

   MAYBE_UNUSED int GRID_DIM = gpuGridSize(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(DEVICE_PTRS,mlst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];
      MAYBE_UNUSED real gxi = 0, gyi = 0, gzi = 0;
      MAYBE_UNUSED real txi = 0, tyi = 0, tzi = 0;

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent private(pgrad)\
                  reduction(+:gxi,gyi,gzi,txi,tyi,tzi)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int offset = (kk + i * n) & (bufsize - 1);
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         real r2 = image2(xr, yr, zr);
         if (r2 <= off2) {
            MAYBE_UNUSED real e;
            pair_mpole<do_e, do_g, EWALD>(r2, xr, yr, zr, 1,          //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, //
               rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
               rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz],
               rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
               rpole[k][mpl_pme_zz], //
               f, aewald, e, pgrad);

            if CONSTEXPR (do_a)
               atomic_add(1, nem, offset);
            if CONSTEXPR (do_e)
               atomic_add(e, em, offset);
            if CONSTEXPR (do_g) {
               gxi += pgrad.frcx;
               gyi += pgrad.frcy;
               gzi += pgrad.frcz;
               atomic_add(-pgrad.frcx, demx, k);
               atomic_add(-pgrad.frcy, demy, k);
               atomic_add(-pgrad.frcz, demz, k);

               txi += pgrad.ttmi[0];
               tyi += pgrad.ttmi[1];
               tzi += pgrad.ttmi[2];
               atomic_add(pgrad.ttmk[0], trqx, k);
               atomic_add(pgrad.ttmk[1], trqy, k);
               atomic_add(pgrad.ttmk[2], trqz, k);

               // virial

               if CONSTEXPR (do_v) {
                  real vxx = -xr * pgrad.frcx;
                  real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
                  real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
                  real vyy = -yr * pgrad.frcy;
                  real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
                  real vzz = -zr * pgrad.frcz;

                  atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_em, offset);
               } // end if (do_v)
            }    // end if (do_g)
         }       // end if (r2 <= off2)
      }          // end for (int kk)

      if CONSTEXPR (do_g) {
         atomic_add(gxi, demx, i);
         atomic_add(gyi, demy, i);
         atomic_add(gzi, demz, i);
         atomic_add(txi, trqx, i);
         atomic_add(tyi, trqy, i);
         atomic_add(tzi, trqz, i);
      }

      // compute the self-energy part of the Ewald summation

      real cii = ci * ci;
      real dii = dix * dix + diy * diy + diz * diz;
      real qii =
         2 * (qixy * qixy + qixz * qixz + qiyz * qiyz) + qixx * qixx + qiyy * qiyy + qizz * qizz;

      if CONSTEXPR (do_e) {
         int offset = i & (bufsize - 1);
         real e = fterm * (cii + aewald_sq_2 * (dii / 3 + 2 * aewald_sq_2 * qii * (real)0.2));
         atomic_add(e, em, offset);
         if CONSTEXPR (do_a)
            atomic_add(1, nem, offset);
      } // end if (do_e)
   }    // end for (int i)

   #pragma acc parallel async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(DEVICE_PTRS,mexclude,mexclude_scale)
   #pragma acc loop independent private(pgrad)
   for (int ii = 0; ii < nmexclude; ++ii) {
      int offset = ii & (bufsize - 1);

      int i = mexclude[ii][0];
      int k = mexclude[ii][1];
      real mscale = mexclude_scale[ii] - 1;

      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];

      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         MAYBE_UNUSED real e;
         pair_mpole<do_e, do_g, NON_EWALD>(                        //
            r2, xr, yr, zr, mscale,                                //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, //
            rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
            rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
            rpole[k][mpl_pme_yz],
            rpole[k][mpl_pme_zz], //
            f, 0, e, pgrad);

         if CONSTEXPR (do_a) {
            if (mscale == -1)
               atomic_add(-1, nem, offset);
         }
         if CONSTEXPR (do_e)
            atomic_add(e, em, offset);
         if CONSTEXPR (do_g) {
            atomic_add(pgrad.frcx, demx, i);
            atomic_add(pgrad.frcy, demy, i);
            atomic_add(pgrad.frcz, demz, i);
            atomic_add(-pgrad.frcx, demx, k);
            atomic_add(-pgrad.frcy, demy, k);
            atomic_add(-pgrad.frcz, demz, k);

            atomic_add(pgrad.ttmi[0], trqx, i);
            atomic_add(pgrad.ttmi[1], trqy, i);
            atomic_add(pgrad.ttmi[2], trqz, i);
            atomic_add(pgrad.ttmk[0], trqx, k);
            atomic_add(pgrad.ttmk[1], trqy, k);
            atomic_add(pgrad.ttmk[2], trqz, k);

            // virial

            if CONSTEXPR (do_v) {
               real vxx = -xr * pgrad.frcx;
               real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               real vyy = -yr * pgrad.frcy;
               real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               real vzz = -zr * pgrad.frcz;

               atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_em, offset);
            } // end if (do_v)
         }    // end if (do_g)
      }
   }
}

void empole_ewald_real_self_acc(int vers)
{
   if (vers == calc::v0)
      empole_ewald_real_self_acc1<calc::V0>();
   else if (vers == calc::v1)
      empole_ewald_real_self_acc1<calc::V1>();
   else if (vers == calc::v3)
      empole_ewald_real_self_acc1<calc::V3>();
   else if (vers == calc::v4)
      empole_ewald_real_self_acc1<calc::V4>();
   else if (vers == calc::v5)
      empole_ewald_real_self_acc1<calc::V5>();
   else if (vers == calc::v6)
      empole_ewald_real_self_acc1<calc::V6>();
}

void empole_ewald_recip_acc(int vers)
{
   if (vers == calc::v0)
      empole_generic_ewald_recip_acc<calc::V0, 0>();
   else if (vers == calc::v1)
      empole_generic_ewald_recip_acc<calc::V1, 0>();
   else if (vers == calc::v3)
      empole_generic_ewald_recip_acc<calc::V3, 0>();
   else if (vers == calc::v4)
      empole_generic_ewald_recip_acc<calc::V4, 0>();
   else if (vers == calc::v5)
      empole_generic_ewald_recip_acc<calc::V5, 0>();
   else if (vers == calc::v6)
      empole_generic_ewald_recip_acc<calc::V6, 0>();
}
}
