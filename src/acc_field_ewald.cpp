#include "acc_add.h"
#include "acc_field_pair.h"
#include "acc_image.h"
#include "couple.h" //
#include "e_mpole.h"
#include "e_polar.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"
#include "pme.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
// see also subroutine udirect1 in induce.f
void dfield_ewald_recip_self(real (*field)[3]) {
  const PMEUnit pu = ppme_unit;
  const real aewald = pu->aewald;
  const real term = REAL_CUBE(aewald) * 4 / 3 / sqrtpi;

  cmp_to_fmp(pu, cmp, fmp);
  grid_mpole(pu, fmp);
  fftfront(pu);
  if (vir_m_handle.valid() && !use_potent(mpole_term))
    pme_conv1(pu, vir_m_handle);
  else
    pme_conv0(pu);
  fftback(pu);
  fphi_mpole(pu, fphi);
  fphi_to_cphi(pu, fphi, cphi);

  #pragma acc parallel loop independent deviceptr(field,cphi,rpole)
  for (int i = 0; i < n; ++i) {
    real dix = rpole[i][mpl_pme_x];
    real diy = rpole[i][mpl_pme_y];
    real diz = rpole[i][mpl_pme_z];
    field[i][0] += (-cphi[i][1] + term * dix);
    field[i][1] += (-cphi[i][2] + term * diy);
    field[i][2] += (-cphi[i][3] + term * diz);
  }
}

// see also subroutine udirect2b / dfield0c in induce.f
void dfield_ewald_real(real (*field)[3], real (*fieldp)[3]) {
  const real off = ewald_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_unit->maxnlst;
  const auto* mlst = mlist_unit.deviceptr();

  const PMEUnit pu = ppme_unit;
  const real aewald = pu->aewald;

#define DFIELD_DPTRS_ x, y, z, box, rpole, thole, pdamp, field, fieldp

  MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
  #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
              deviceptr(DFIELD_DPTRS_,mlst)
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
    real pdi = pdamp[i];
    real pti = thole[i];
    real gxi = 0, gyi = 0, gzi = 0;
    real txi = 0, tyi = 0, tzi = 0;

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop vector independent reduction(+:gxi,gyi,gzi,txi,tyi,tzi)
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        FieldPair pairf;
        dfield_pair_acc<elec_t::ewald>(         //
            r2, xr, yr, zr,                     //
            1, 1, aewald,                       //
            ci, dix, diy, diz,                  //
            qixx, qixy, qixz, qiyy, qiyz, qizz, //
            pdi, pti,                           //
            rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
            rpole[k][mpl_pme_z], //
            rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz],
            rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], //
            pdamp[k], thole[k],                                               //
            pairf);

        gxi += pairf.fid[0];
        gyi += pairf.fid[1];
        gzi += pairf.fid[2];
        txi += pairf.fip[0];
        tyi += pairf.fip[1];
        tzi += pairf.fip[2];

        atomic_add_value(pairf.fkd[0], &field[k][0]);
        atomic_add_value(pairf.fkd[1], &field[k][1]);
        atomic_add_value(pairf.fkd[2], &field[k][2]);
        atomic_add_value(pairf.fkp[0], &fieldp[k][0]);
        atomic_add_value(pairf.fkp[1], &fieldp[k][1]);
        atomic_add_value(pairf.fkp[2], &fieldp[k][2]);
      }
    } // end for (int kk)

    atomic_add_value(gxi, &field[i][0]);
    atomic_add_value(gyi, &field[i][1]);
    atomic_add_value(gzi, &field[i][2]);
    atomic_add_value(txi, &fieldp[i][0]);
    atomic_add_value(tyi, &fieldp[i][1]);
    atomic_add_value(tzi, &fieldp[i][2]);
  } // end for (int i)

  #pragma acc parallel deviceptr(DFIELD_DPTRS_,dpexclude_,dpexclude_scale_)
  #pragma acc loop independent
  for (int ii = 0; ii < ndpexclude_; ++ii) {
    int i = dpexclude_[ii][0];
    int k = dpexclude_[ii][1];
    real dscale = dpexclude_scale_[ii][0];
    real pscale = dpexclude_scale_[ii][1];

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
    real pdi = pdamp[i];
    real pti = thole[i];

    real xr = x[k] - xi;
    real yr = y[k] - yi;
    real zr = z[k] - zi;

    image(xr, yr, zr, box);
    real r2 = xr * xr + yr * yr + zr * zr;

    FieldPair pairf;
    dfield_pair_acc<elec_t::coulomb>(       //
        r2, xr, yr, zr,                     //
        dscale, pscale, 0,                  //
        ci, dix, diy, diz,                  //
        qixx, qixy, qixz, qiyy, qiyz, qizz, //
        pdi, pti,                           //
        rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
        rpole[k][mpl_pme_z], //
        rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz],
        rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], //
        pdamp[k], thole[k],                                               //
        pairf);

    atomic_add_value(pairf.fid[0], &field[i][0]);
    atomic_add_value(pairf.fid[1], &field[i][1]);
    atomic_add_value(pairf.fid[2], &field[i][2]);
    atomic_add_value(pairf.fip[0], &fieldp[i][0]);
    atomic_add_value(pairf.fip[1], &fieldp[i][1]);
    atomic_add_value(pairf.fip[2], &fieldp[i][2]);

    atomic_add_value(pairf.fkd[0], &field[k][0]);
    atomic_add_value(pairf.fkd[1], &field[k][1]);
    atomic_add_value(pairf.fkd[2], &field[k][2]);
    atomic_add_value(pairf.fkp[0], &fieldp[k][0]);
    atomic_add_value(pairf.fkp[1], &fieldp[k][1]);
    atomic_add_value(pairf.fkp[2], &fieldp[k][2]);
  }
}

void dfield_ewald(real (*field)[3], real (*fieldp)[3]) {
  device_array::zero(n, field, fieldp);

  dfield_ewald_recip_self(field);
  device_array::copy(n, fieldp, field);

  dfield_ewald_real(field, fieldp);
}

// see also subroutine umutual1 in induce.f
void ufield_ewald_recip_self(const real (*uind)[3], const real (*uinp)[3],
                             real (*field)[3], real (*fieldp)[3]) {
  const PMEUnit pu = ppme_unit;
  const auto& st = *pu;
  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const real aewald = st.aewald;

  cuind_to_fuind(pu, uind, uinp, fuind, fuinp);
  grid_uind(pu, fuind, fuinp);
  fftfront(pu);
  // TODO: store vs. recompute qfac
  pme_conv0(pu);
  fftback(pu);
  fphi_uind2(pu, fdip_phi1, fdip_phi2);

  real a[3][3];

  const real term = REAL_CUBE(aewald) * 4 / 3 / sqrtpi;

  #pragma acc parallel loop independent deviceptr(box,\
              field,fieldp,uind,uinp,\
              fdip_phi1,fdip_phi2)\
              private(a[0:3][0:3])
  for (int i = 0; i < n; ++i) {
    a[0][0] = nfft1 * box->recip[0][0];
    a[1][0] = nfft2 * box->recip[1][0];
    a[2][0] = nfft3 * box->recip[2][0];
    a[0][1] = nfft1 * box->recip[0][1];
    a[1][1] = nfft2 * box->recip[1][1];
    a[2][1] = nfft3 * box->recip[2][1];
    a[0][2] = nfft1 * box->recip[0][2];
    a[1][2] = nfft2 * box->recip[1][2];
    a[2][2] = nfft3 * box->recip[2][2];

    #pragma acc loop independent
    for (int j = 0; j < 3; ++j) {
      real df1 = a[0][j] * fdip_phi1[i][1] + a[1][j] * fdip_phi1[i][2] +
          a[2][j] * fdip_phi1[i][3];
      real df2 = a[0][j] * fdip_phi2[i][1] + a[1][j] * fdip_phi2[i][2] +
          a[2][j] * fdip_phi2[i][3];
      field[i][j] += (term * uind[i][j] - df1);
      fieldp[i][j] += (term * uinp[i][j] - df2);
    }
  }
}

void ufield_ewald_real(const real (*uind)[3], const real (*uinp)[3],
                       real (*field)[3], real (*fieldp)[3]) {
  const real off = ewald_switch_cut;
  const real off2 = off * off;
  const int maxnlst = mlist_unit->maxnlst;
  const auto* mlst = mlist_unit.deviceptr();

  const auto* polargroup = polargroup_unit.deviceptr();

  auto bufsize = EnergyBuffer::calc_size(n);

  static std::vector<real> uscalebuf;
  uscalebuf.resize(n, 1);
  real* uscale = uscalebuf.data();

  const PMEUnit pu = ppme_unit;
  const real aewald = pu->aewald;
  const real aesq2 = 2 * aewald * aewald;
  const real aesq2n = (aewald > 0 ? 1 / (sqrtpi * aewald) : 0);

  real bn[3];

  #pragma acc parallel num_gangs(bufsize)\
              deviceptr(x,y,z,box,polargroup,mlst,\
              thole,pdamp,uind,uinp,field,fieldp)\
              firstprivate(uscale[0:n])
  #pragma acc loop gang independent
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int np11i = polargroup->np11[i];
    const int np12i = polargroup->np12[i];
    const int np13i = polargroup->np13[i];
    const int np14i = polargroup->np14[i];

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      uscale[polargroup->ip11[i][j]] = u1scale;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      uscale[polargroup->ip12[i][j]] = u2scale;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      uscale[polargroup->ip13[i][j]] = u3scale;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      uscale[polargroup->ip14[i][j]] = u4scale;

    real xi = x[i];
    real yi = y[i];
    real zi = z[i];

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent private(bn[0:3])
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real pdi = pdamp[i];
      real pti = thole[i];

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        real r = REAL_SQRT(r2);
        real rr1 = REAL_RECIP(r);
        real rr2 = rr1 * rr1;
        real rr3 = rr1 * rr2;
        real rr5 = 3 * rr3 * rr2;

        real ralpha = aewald * r;
        bn[0] = REAL_ERFC(ralpha) * rr1;
        real exp2a = REAL_EXP(-REAL_SQ(ralpha));
        real aefac = aesq2n;
        #pragma acc loop seq
        for (int j = 1; j <= 2; ++j) {
          aefac *= aesq2;
          bn[j] = ((j + j - 1) * bn[j - 1] + aefac * exp2a) * rr2;
        }

        // if use_thole
        real scale3 = 1;
        real scale5 = 1;
        real damp = pdi * pdamp[k];
        if (damp != 0) {
          real pgamma = REAL_MIN(pti, thole[k]);
          damp = -pgamma * REAL_CUBE(r / damp);
          if (damp > -50) {
            real expdamp = REAL_EXP(damp);
            scale3 = 1 - expdamp;
            scale5 = 1 - expdamp * (1 - damp);
          }
        }
        real scalek, bcn1, bcn2;

        scalek = uscale[k];
        bcn1 = bn[1] - (1 - scalek * scale3) * rr3;
        bcn2 = bn[2] - (1 - scalek * scale5) * rr5;

        real dlocal1 = bcn2 * xr * xr - bcn1;
        real dlocal2 = bcn2 * xr * yr;
        real dlocal3 = bcn2 * xr * zr;
        real dlocal4 = bcn2 * yr * yr - bcn1;
        real dlocal5 = bcn2 * yr * zr;
        real dlocal6 = bcn2 * zr * zr - bcn1;

        real fid1 =
            dlocal1 * uind[k][0] + dlocal2 * uind[k][1] + dlocal3 * uind[k][2];
        real fid2 =
            dlocal2 * uind[k][0] + dlocal4 * uind[k][1] + dlocal5 * uind[k][2];
        real fid3 =
            dlocal3 * uind[k][0] + dlocal5 * uind[k][1] + dlocal6 * uind[k][2];
        real fkd1 =
            dlocal1 * uind[i][0] + dlocal2 * uind[i][1] + dlocal3 * uind[i][2];
        real fkd2 =
            dlocal2 * uind[i][0] + dlocal4 * uind[i][1] + dlocal5 * uind[i][2];
        real fkd3 =
            dlocal3 * uind[i][0] + dlocal5 * uind[i][1] + dlocal6 * uind[i][2];

        real fip1 =
            dlocal1 * uinp[k][0] + dlocal2 * uinp[k][1] + dlocal3 * uinp[k][2];
        real fip2 =
            dlocal2 * uinp[k][0] + dlocal4 * uinp[k][1] + dlocal5 * uinp[k][2];
        real fip3 =
            dlocal3 * uinp[k][0] + dlocal5 * uinp[k][1] + dlocal6 * uinp[k][2];
        real fkp1 =
            dlocal1 * uinp[i][0] + dlocal2 * uinp[i][1] + dlocal3 * uinp[i][2];
        real fkp2 =
            dlocal2 * uinp[i][0] + dlocal4 * uinp[i][1] + dlocal5 * uinp[i][2];
        real fkp3 =
            dlocal3 * uinp[i][0] + dlocal5 * uinp[i][1] + dlocal6 * uinp[i][2];

        #pragma acc atomic update
        field[i][0] += fid1;
        #pragma acc atomic update
        field[i][1] += fid2;
        #pragma acc atomic update
        field[i][2] += fid3;
        #pragma acc atomic update
        field[k][0] += fkd1;
        #pragma acc atomic update
        field[k][1] += fkd2;
        #pragma acc atomic update
        field[k][2] += fkd3;

        #pragma acc atomic update
        fieldp[i][0] += fip1;
        #pragma acc atomic update
        fieldp[i][1] += fip2;
        #pragma acc atomic update
        fieldp[i][2] += fip3;
        #pragma acc atomic update
        fieldp[k][0] += fkp1;
        #pragma acc atomic update
        fieldp[k][1] += fkp2;
        #pragma acc atomic update
        fieldp[k][2] += fkp3;
        // end if use_thole
      }
    } // end for (int kk)

    // reset exclusion coefficients for connected atoms

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      uscale[polargroup->ip11[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      uscale[polargroup->ip12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      uscale[polargroup->ip13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      uscale[polargroup->ip14[i][j]] = 1;
  } // end for (int i)
}

void ufield_ewald(const real (*uind)[3], const real (*uinp)[3],
                  real (*field)[3], real (*fieldp)[3]) {
  device_array::zero(n, field, fieldp);

  ufield_ewald_recip_self(uind, uinp, field, fieldp);
  ufield_ewald_real(uind, uinp, field, fieldp);
}
TINKER_NAMESPACE_END
