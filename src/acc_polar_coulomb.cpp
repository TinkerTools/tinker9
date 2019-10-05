#include "acc_common.h"
#include "acc_polar_pair.h"
#include "e_polar.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void epolar_coulomb_tmpl(const real (*uind)[3], const real (*uinp)[3]) {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  if_constexpr(do_g) device_array::zero(n, ufld, dufld);

  if_constexpr(do_e && !do_a) epolar0_dotprod(uind, udirp);
  static_assert(do_g || do_a,
                "Do not use this template for the energy-only version.");

  const real off = switch_off(switch_mpole);
  const real off2 = off * off;
  const int maxnlst = mlist_unit->maxnlst;
  const auto* mlst = mlist_unit.deviceptr();

  auto* nep = ep_handle.ne()->buffer();
  auto* ep = ep_handle.e()->buffer();
  auto* vir_ep = ep_handle.vir()->buffer();
  auto bufsize = ep_handle.buffer_size();

  const real f = 0.5 * electric / dielec;

#define POLAR_DPTRS_                                                           \
  x, y, z, gx, gy, gz, box, rpole, thole, pdamp, uind, uinp, nep, ep, vir_ep,  \
      ufld, dufld

  MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
  #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
              deviceptr(POLAR_DPTRS_,mlst)
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
    real uix = uind[i][0];
    real uiy = uind[i][1];
    real uiz = uind[i][2];
    real pdi = pdamp[i];
    real pti = thole[i];
    MAYBE_UNUSED real uixp, uiyp, uizp;
    if_constexpr(do_g) {
      uixp = uinp[i][0];
      uiyp = uinp[i][1];
      uizp = uinp[i][2];
    }
    real gxi = 0, gyi = 0, gzi = 0;
    real txi = 0, tyi = 0, tzi = 0;
    real du0 = 0, du1 = 0, du2 = 0, du3 = 0, du4 = 0, du5 = 0;

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop vector independent\
                reduction(+:gxi,gyi,gzi,txi,tyi,tzi,du0,du1,du2,du3,du4,du5)
    for (int kk = 0; kk < nmlsti; ++kk) {
      int offset = kk & (bufsize - 1);
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        real ck = rpole[k][mpl_pme_0];
        real dkx = rpole[k][mpl_pme_x];
        real dky = rpole[k][mpl_pme_y];
        real dkz = rpole[k][mpl_pme_z];
        real qkxx = rpole[k][mpl_pme_xx];
        real qkxy = rpole[k][mpl_pme_xy];
        real qkxz = rpole[k][mpl_pme_xz];
        real qkyy = rpole[k][mpl_pme_yy];
        real qkyz = rpole[k][mpl_pme_yz];
        real qkzz = rpole[k][mpl_pme_zz];
        real ukx = uind[k][0];
        real uky = uind[k][1];
        real ukz = uind[k][2];
        MAYBE_UNUSED real ukxp, ukyp, ukzp;
        if_constexpr(do_g) {
          ukxp = uinp[k][0];
          ukyp = uinp[k][1];
          ukzp = uinp[k][2];
        }

        MAYBE_UNUSED real e;
        MAYBE_UNUSED PolarPairGrad pgrad;
        epolar_pair_acc<USE, elec_t::coulomb>(  //
            r2, xr, yr, zr,                     //
            f, 1, 1, 1, 0,                      //
            ci, dix, diy, diz,                  //
            qixx, qixy, qixz, qiyy, qiyz, qizz, //
            //
            uix, uiy, uiz,    //
            uixp, uiyp, uizp, //
            pdi, pti,         //
            //
            ck, dkx, dky, dkz,                  //
            qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, //
            ukx, uky, ukz,                      //
            ukxp, ukyp, ukzp,                   //
            pdamp[k], thole[k],                 //
            e, pgrad);

        if_constexpr(do_a && do_e) {
          atomic_add_value(1, nep, offset);
          atomic_add_value(e, ep, offset);
        }

        if_constexpr(do_g) {
          gxi += pgrad.frcx;
          gyi += pgrad.frcy;
          gzi += pgrad.frcz;
          atomic_add_value(-pgrad.frcx, gx, k);
          atomic_add_value(-pgrad.frcy, gy, k);
          atomic_add_value(-pgrad.frcz, gz, k);

          txi += pgrad.ufldi[0];
          tyi += pgrad.ufldi[1];
          tzi += pgrad.ufldi[2];
          atomic_add_value(pgrad.ufldk[0], &ufld[k][0]);
          atomic_add_value(pgrad.ufldk[1], &ufld[k][1]);
          atomic_add_value(pgrad.ufldk[2], &ufld[k][2]);

          du0 += pgrad.dufldi[0];
          du1 += pgrad.dufldi[1];
          du2 += pgrad.dufldi[2];
          du3 += pgrad.dufldi[3];
          du4 += pgrad.dufldi[4];
          du5 += pgrad.dufldi[5];
          atomic_add_value(pgrad.dufldk[0], &dufld[k][0]);
          atomic_add_value(pgrad.dufldk[1], &dufld[k][1]);
          atomic_add_value(pgrad.dufldk[2], &dufld[k][2]);
          atomic_add_value(pgrad.dufldk[3], &dufld[k][3]);
          atomic_add_value(pgrad.dufldk[4], &dufld[k][4]);
          atomic_add_value(pgrad.dufldk[5], &dufld[k][5]);

          if_constexpr(do_v) {
            real vxx = -xr * pgrad.frcx;
            real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
            real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
            real vyy = -yr * pgrad.frcy;
            real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
            real vzz = -zr * pgrad.frcz;

            atomic_add_value(vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, offset);
          }
        }
      }
    } // end for (int kk)

    if_constexpr(do_g) {
      atomic_add_value(gxi, gx, i);
      atomic_add_value(gyi, gy, i);
      atomic_add_value(gzi, gz, i);
      atomic_add_value(txi, &ufld[i][0]);
      atomic_add_value(tyi, &ufld[i][1]);
      atomic_add_value(tzi, &ufld[i][2]);
      atomic_add_value(du0, &dufld[i][0]);
      atomic_add_value(du1, &dufld[i][1]);
      atomic_add_value(du2, &dufld[i][2]);
      atomic_add_value(du3, &dufld[i][3]);
      atomic_add_value(du4, &dufld[i][4]);
      atomic_add_value(du5, &dufld[i][5]);
    }
  } // end for (int i)

  #pragma acc parallel deviceptr(POLAR_DPTRS_,dpuexclude_,dpuexclude_scale_)
  #pragma acc loop independent
  for (int ii = 0; ii < ndpuexclude_; ++ii) {
    int offset = ii & (bufsize - 1);

    int i = dpuexclude_[ii][0];
    int k = dpuexclude_[ii][1];
    real dscale = dpuexclude_scale_[ii][0];
    real pscale = dpuexclude_scale_[ii][1];
    real uscale = dpuexclude_scale_[ii][2];

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
    real uix = uind[i][0];
    real uiy = uind[i][1];
    real uiz = uind[i][2];
    real pdi = pdamp[i];
    real pti = thole[i];
    MAYBE_UNUSED real uixp, uiyp, uizp;
    if_constexpr(do_g) {
      uixp = uinp[i][0];
      uiyp = uinp[i][1];
      uizp = uinp[i][2];
    }

    real xr = x[k] - xi;
    real yr = y[k] - yi;
    real zr = z[k] - zi;

    image(xr, yr, zr, box);
    real r2 = xr * xr + yr * yr + zr * zr;
    if (r2 <= off2) {
      real ck = rpole[k][mpl_pme_0];
      real dkx = rpole[k][mpl_pme_x];
      real dky = rpole[k][mpl_pme_y];
      real dkz = rpole[k][mpl_pme_z];
      real qkxx = rpole[k][mpl_pme_xx];
      real qkxy = rpole[k][mpl_pme_xy];
      real qkxz = rpole[k][mpl_pme_xz];
      real qkyy = rpole[k][mpl_pme_yy];
      real qkyz = rpole[k][mpl_pme_yz];
      real qkzz = rpole[k][mpl_pme_zz];
      real ukx = uind[k][0];
      real uky = uind[k][1];
      real ukz = uind[k][2];
      MAYBE_UNUSED real ukxp, ukyp, ukzp;
      if_constexpr(do_g) {
        ukxp = uinp[k][0];
        ukyp = uinp[k][1];
        ukzp = uinp[k][2];
      }

      MAYBE_UNUSED real e;
      MAYBE_UNUSED PolarPairGrad pgrad;
      epolar_pair_acc<USE, elec_t::coulomb>(  //
          r2, xr, yr, zr,                     //
          f, dscale, pscale, uscale, 0,       //
          ci, dix, diy, diz,                  //
          qixx, qixy, qixz, qiyy, qiyz, qizz, //
          //
          uix, uiy, uiz,    //
          uixp, uiyp, uizp, //
          pdi, pti,         //
          //
          ck, dkx, dky, dkz,                  //
          qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, //
          ukx, uky, ukz,                      //
          ukxp, ukyp, ukzp,                   //
          pdamp[k], thole[k],                 //
          e, pgrad);

      if_constexpr(do_a && do_e) {
        if (pscale == -1)
          atomic_add_value(-1, nep, offset);
        atomic_add_value(e, ep, offset);
      }

      if_constexpr(do_g) {
        atomic_add_value(pgrad.frcx, gx, i);
        atomic_add_value(pgrad.frcy, gy, i);
        atomic_add_value(pgrad.frcz, gz, i);
        atomic_add_value(-pgrad.frcx, gx, k);
        atomic_add_value(-pgrad.frcy, gy, k);
        atomic_add_value(-pgrad.frcz, gz, k);

        atomic_add_value(pgrad.ufldi[0], &ufld[i][0]);
        atomic_add_value(pgrad.ufldi[1], &ufld[i][1]);
        atomic_add_value(pgrad.ufldi[2], &ufld[i][2]);
        atomic_add_value(pgrad.ufldk[0], &ufld[k][0]);
        atomic_add_value(pgrad.ufldk[1], &ufld[k][1]);
        atomic_add_value(pgrad.ufldk[2], &ufld[k][2]);

        atomic_add_value(pgrad.dufldi[0], &dufld[i][0]);
        atomic_add_value(pgrad.dufldi[1], &dufld[i][1]);
        atomic_add_value(pgrad.dufldi[2], &dufld[i][2]);
        atomic_add_value(pgrad.dufldi[3], &dufld[i][3]);
        atomic_add_value(pgrad.dufldi[4], &dufld[i][4]);
        atomic_add_value(pgrad.dufldi[5], &dufld[i][5]);
        atomic_add_value(pgrad.dufldk[0], &dufld[k][0]);
        atomic_add_value(pgrad.dufldk[1], &dufld[k][1]);
        atomic_add_value(pgrad.dufldk[2], &dufld[k][2]);
        atomic_add_value(pgrad.dufldk[3], &dufld[k][3]);
        atomic_add_value(pgrad.dufldk[4], &dufld[k][4]);
        atomic_add_value(pgrad.dufldk[5], &dufld[k][5]);

        if_constexpr(do_v) {
          real vxx = -xr * pgrad.frcx;
          real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
          real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
          real vyy = -yr * pgrad.frcy;
          real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
          real vzz = -zr * pgrad.frcz;

          atomic_add_value(vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, offset);
        }
      }
    }
  }

  // torque

  if_constexpr(do_g) {
    #pragma acc parallel loop independent\
                deviceptr(rpole,trqx,trqy,trqz,ufld,dufld)
    for (int i = 0; i < n; ++i) {
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];

      real tep1 = diz * ufld[i][1] - diy * ufld[i][2] + qixz * dufld[i][1] -
          qixy * dufld[i][3] + 2 * qiyz * (dufld[i][2] - dufld[i][5]) +
          (qizz - qiyy) * dufld[i][4];
      real tep2 = dix * ufld[i][2] - diz * ufld[i][0] - qiyz * dufld[i][1] +
          qixy * dufld[i][4] + 2 * qixz * (dufld[i][5] - dufld[i][0]) +
          (qixx - qizz) * dufld[i][3];
      real tep3 = diy * ufld[i][0] - dix * ufld[i][1] + qiyz * dufld[i][3] -
          qixz * dufld[i][4] + 2 * qixy * (dufld[i][0] - dufld[i][2]) +
          (qiyy - qixx) * dufld[i][1];

      trqx[i] += tep1;
      trqy[i] += tep2;
      trqz[i] += tep3;
    }
  }
}

void epolar_coulomb(int vers) {
  if (vers == calc::v0) {
    induce(uind, uinp);
    epolar0_dotprod(uind, udirp);
  } else if (vers == calc::v1) {
    induce(uind, uinp);
    epolar_coulomb_tmpl<calc::v1>(uind, uinp);
  } else if (vers == calc::v3) {
    induce(uind, uinp);
    epolar_coulomb_tmpl<calc::v3>(uind, uinp);
  } else if (vers == calc::v4) {
    induce(uind, uinp);
    epolar_coulomb_tmpl<calc::v4>(uind, uinp);
  } else if (vers == calc::v5) {
    induce(uind, uinp);
    epolar_coulomb_tmpl<calc::v5>(uind, uinp);
  } else if (vers == calc::v6) {
    induce(uind, uinp);
    epolar_coulomb_tmpl<calc::v6>(uind, uinp);
  }
}
TINKER_NAMESPACE_END
