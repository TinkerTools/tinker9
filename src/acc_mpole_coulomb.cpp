#include "acc_add.h"
#include "e_mpole.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"
#include "seq_image.h"
#include "seq_pair_mpole.h"
#include "seq_switch.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void empole_coulomb_tmpl()
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   static_assert(do_v ? do_g : true, "");
   static_assert(do_a ? do_e : true, "");

   const real f = electric / dielec;

   const real off = switch_off(switch_mpole);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   auto bufsize = buffer_size();

#define DEVICE_PTRS_                                                           \
   x, y, z, gx, gy, gz, box, rpole, nem, em, vir_em, trqx, trqy, trqz

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(DEVICE_PTRS_,mlst)
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
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi,txi,tyi,tzi)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int offset = (kk + i * n) & (bufsize - 1);
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (r2 <= off2) {
            MAYBE_UNUSED real e;
            MAYBE_UNUSED PairMPoleGrad pgrad;
            pair_mpole<USE, elec_t::coulomb>(                         //
               r2, xr, yr, zr, 1,                                     //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, //
               rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
               rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy],
               rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
               rpole[k][mpl_pme_zz], //
               f, 0, e, pgrad);

            if_constexpr(do_a) atomic_add_value(1, nem, offset);
            if_constexpr(do_e) atomic_add_value(e, em, offset);
            if_constexpr(do_g)
            {
               gxi += pgrad.frcx;
               gyi += pgrad.frcy;
               gzi += pgrad.frcz;
               atomic_add_value(-pgrad.frcx, gx, k);
               atomic_add_value(-pgrad.frcy, gy, k);
               atomic_add_value(-pgrad.frcz, gz, k);

               txi += pgrad.ttmi[0];
               tyi += pgrad.ttmi[1];
               tzi += pgrad.ttmi[2];
               atomic_add_value(pgrad.ttmk[0], trqx, k);
               atomic_add_value(pgrad.ttmk[1], trqy, k);
               atomic_add_value(pgrad.ttmk[2], trqz, k);

               // virial

               if_constexpr(do_v)
               {
                  real vxx = -xr * pgrad.frcx;
                  real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
                  real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
                  real vyy = -yr * pgrad.frcy;
                  real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
                  real vzz = -zr * pgrad.frcz;

                  atomic_add_value(vxx, vxy, vxz, vyy, vyz, vzz, vir_em,
                                   offset);
               } // end if (do_v)
            }    // end if (do_g)
         }       // end if (r2 <= off2)
      }          // end for (int kk)

      if_constexpr(do_g)
      {
         atomic_add_value(gxi, gx, i);
         atomic_add_value(gyi, gy, i);
         atomic_add_value(gzi, gz, i);
         atomic_add_value(txi, trqx, i);
         atomic_add_value(tyi, trqy, i);
         atomic_add_value(tzi, trqz, i);
      }
   } // end for (int i)

   #pragma acc parallel deviceptr(DEVICE_PTRS_,mexclude_,mexclude_scale_)
   #pragma acc loop independent
   for (int ii = 0; ii < nmexclude_; ++ii) {
      int offset = ii & (bufsize - 1);

      int i = mexclude_[ii][0];
      int k = mexclude_[ii][1];
      real mscale = mexclude_scale_[ii];

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

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
         MAYBE_UNUSED real e;
         MAYBE_UNUSED PairMPoleGrad pgrad;
         pair_mpole<USE, elec_t::coulomb>(                         //
            r2, xr, yr, zr, mscale,                                //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, //
            rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
            rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy],
            rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
            rpole[k][mpl_pme_zz], //
            f, 0, e, pgrad);

         if_constexpr(do_a) if (mscale == -1) atomic_add_value(-1, nem, offset);
         if_constexpr(do_e) atomic_add_value(e, em, offset);
         if_constexpr(do_g)
         {
            atomic_add_value(pgrad.frcx, gx, i);
            atomic_add_value(pgrad.frcy, gy, i);
            atomic_add_value(pgrad.frcz, gz, i);
            atomic_add_value(-pgrad.frcx, gx, k);
            atomic_add_value(-pgrad.frcy, gy, k);
            atomic_add_value(-pgrad.frcz, gz, k);

            atomic_add_value(pgrad.ttmi[0], trqx, i);
            atomic_add_value(pgrad.ttmi[1], trqy, i);
            atomic_add_value(pgrad.ttmi[2], trqz, i);
            atomic_add_value(pgrad.ttmk[0], trqx, k);
            atomic_add_value(pgrad.ttmk[1], trqy, k);
            atomic_add_value(pgrad.ttmk[2], trqz, k);

            // virial

            if_constexpr(do_v)
            {
               real vxx = -xr * pgrad.frcx;
               real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               real vyy = -yr * pgrad.frcy;
               real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               real vzz = -zr * pgrad.frcz;

               atomic_add_value(vxx, vxy, vxz, vyy, vyz, vzz, vir_em, offset);
            } // end if (do_v)
         }    // end if (do_g)
      }
   }
}

void empole_coulomb(int vers)
{
   if (vers == calc::v0)
      empole_coulomb_tmpl<calc::v0>();
   else if (vers == calc::v1)
      empole_coulomb_tmpl<calc::v1>();
   else if (vers == calc::v3)
      empole_coulomb_tmpl<calc::v3>();
   else if (vers == calc::v4)
      empole_coulomb_tmpl<calc::v4>();
   else if (vers == calc::v5)
      empole_coulomb_tmpl<calc::v5>();
   else if (vers == calc::v6)
      empole_coulomb_tmpl<calc::v6>();
}
TINKER_NAMESPACE_END
