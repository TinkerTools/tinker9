#include "acc_seq.h"
#include "array.h"
#include "couple.h"
#include "e_polar.h"
#include "md.h"
#include "nblist.h"

TINKER_NAMESPACE_BEGIN
// see also subroutine dfield0b in induce.f
void dfield_coulomb(real* gpu_field, real* gpu_fieldp) {
  zero_array(gpu_field, 3 * n);
  zero_array(gpu_fieldp, 3 * n);

  real(*field)[3] = reinterpret_cast<real(*)[3]>(gpu_field);
  real(*fieldp)[3] = reinterpret_cast<real(*)[3]>(gpu_fieldp);

  const real off = mpole_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_unit->maxnlst;
  const auto* mlst = mlist_unit.deviceptr();

  const auto* coupl = couple_unit.deviceptr();
  const auto* polargroup = polargroup_unit.deviceptr();

  static std::vector<real> pscalebuf;
  static std::vector<real> dscalebuf;
  pscalebuf.resize(n, 1);
  dscalebuf.resize(n, 1);
  real* pscale = pscalebuf.data();
  real* dscale = dscalebuf.data();

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,box,coupl,polargroup,mlst,\
              rpole,thole,pdamp,\
              field,fieldp)\
              firstprivate(pscale[0:n],dscale[0:n])
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int n12i = coupl->n12[i];
    const int n13i = coupl->n13[i];
    const int n14i = coupl->n14[i];
    const int n15i = coupl->n15[i];

    const int np11i = polargroup->np11[i];
    const int np12i = polargroup->np12[i];
    const int np13i = polargroup->np13[i];
    const int np14i = polargroup->np14[i];

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j) {
      int ij = coupl->i12[i][j];
      pscale[ij] = p2scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p2iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j) {
      int ij = coupl->i13[i][j];
      pscale[ij] = p3scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p3iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j) {
      int ij = coupl->i14[i][j];
      pscale[ij] = p4scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p4iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j) {
      int ij = coupl->i15[i][j];
      pscale[ij] = p5scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p5iscale;
    }

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      dscale[polargroup->ip11[i][j]] = d1scale;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      dscale[polargroup->ip12[i][j]] = d2scale;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      dscale[polargroup->ip13[i][j]] = d3scale;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      dscale[polargroup->ip14[i][j]] = d4scale;

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

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent
    for (int kk = 0; kk < nmlsti; ++kk) {
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

        real dir = dix * xr + diy * yr + diz * zr;
        real qix = qixx * xr + qixy * yr + qixz * zr;
        real qiy = qixy * xr + qiyy * yr + qiyz * zr;
        real qiz = qixz * xr + qiyz * yr + qizz * zr;
        real qir = qix * xr + qiy * yr + qiz * zr;
        real dkr = dkx * xr + dky * yr + dkz * zr;
        real qkx = qkxx * xr + qkxy * yr + qkxz * zr;
        real qky = qkxy * xr + qkyy * yr + qkyz * zr;
        real qkz = qkxz * xr + qkyz * yr + qkzz * zr;
        real qkr = qkx * xr + qky * yr + qkz * zr;

        real r = REAL_SQRT(r2);

        // if use_thole
        real scale3 = 1;
        real scale5 = 1;
        real scale7 = 1;
        real damp = pdi * pdamp[k];
        if (damp != 0) {
          real pgamma = REAL_MIN(pti, thole[k]);
          damp = -pgamma * REAL_CUBE(r / damp);
          if (damp > -50) {
            real expdamp = REAL_EXP(damp);
            scale3 = 1 - expdamp;
            scale5 = 1 - expdamp * (1 - damp);
            scale7 = 1 - expdamp * (1 - damp + (real)0.6 * REAL_SQ(damp));
          }
        }

        real rr1 = REAL_RECIP(r);
        real rr2 = rr1 * rr1;

        real rr3 = scale3 * rr1 * rr2;
        real rr5 = 3 * scale5 * rr1 * rr2 * rr2;
        real rr7 = 15 * scale7 * rr1 * rr2 * rr2 * rr2;

        real fid1, fid2, fid3;
        real fkd1, fkd2, fkd3;
        fid1 = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkx +
            2 * rr5 * qkx;
        fid2 = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dky +
            2 * rr5 * qky;
        fid3 = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkz +
            2 * rr5 * qkz;
        fkd1 =
            xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * dix - 2 * rr5 * qix;
        fkd2 =
            yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diy - 2 * rr5 * qiy;
        fkd3 =
            zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diz - 2 * rr5 * qiz;
        // end if use_thole

        real dscalek = dscale[k];
        #pragma acc atomic update
        field[i][0] += fid1 * dscalek;
        #pragma acc atomic update
        field[i][1] += fid2 * dscalek;
        #pragma acc atomic update
        field[i][2] += fid3 * dscalek;
        #pragma acc atomic update
        field[k][0] += fkd1 * dscalek;
        #pragma acc atomic update
        field[k][1] += fkd2 * dscalek;
        #pragma acc atomic update
        field[k][2] += fkd3 * dscalek;

        real pscalek = pscale[k];
        #pragma acc atomic update
        fieldp[i][0] += fid1 * pscalek;
        #pragma acc atomic update
        fieldp[i][1] += fid2 * pscalek;
        #pragma acc atomic update
        fieldp[i][2] += fid3 * pscalek;
        #pragma acc atomic update
        fieldp[k][0] += fkd1 * pscalek;
        #pragma acc atomic update
        fieldp[k][1] += fkd2 * pscalek;
        #pragma acc atomic update
        fieldp[k][2] += fkd3 * pscalek;
      }
    } // end for (int kk)

    // reset exclusion coefficients for connected atoms

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      pscale[coupl->i12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      pscale[coupl->i13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      pscale[coupl->i14[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      pscale[coupl->i15[i][j]] = 1;

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      dscale[polargroup->ip11[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      dscale[polargroup->ip12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      dscale[polargroup->ip13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      dscale[polargroup->ip14[i][j]] = 1;
  } // end for (int i)
}

// see also subroutine ufield0b in induce.f
void ufield_coulomb(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
                    real* gpu_fieldp) {
  zero_array(gpu_field, 3 * n);
  zero_array(gpu_fieldp, 3 * n);

  const real(*uind)[3] = reinterpret_cast<const real(*)[3]>(gpu_uind);
  const real(*uinp)[3] = reinterpret_cast<const real(*)[3]>(gpu_uinp);
  real(*field)[3] = reinterpret_cast<real(*)[3]>(gpu_field);
  real(*fieldp)[3] = reinterpret_cast<real(*)[3]>(gpu_fieldp);

  const real off = mpole_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_unit->maxnlst;
  const auto* mlst = mlist_unit.deviceptr();

  const auto* polargroup = polargroup_unit.deviceptr();

  static std::vector<real> uscalebuf;
  uscalebuf.resize(n, 1);
  real* uscale = uscalebuf.data();

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,box,polargroup,mlst,\
              thole,pdamp,uind,uinp,field,fieldp)\
              firstprivate(uscale[0:n])
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
    #pragma acc loop independent
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real dix = uind[i][0];
      real diy = uind[i][1];
      real diz = uind[i][2];
      real pix = uinp[i][0];
      real piy = uinp[i][1];
      real piz = uinp[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        real dkx = uind[k][0];
        real dky = uind[k][1];
        real dkz = uind[k][2];
        real pkx = uinp[k][0];
        real pky = uinp[k][1];
        real pkz = uinp[k][2];

        real dir = dix * xr + diy * yr + diz * zr;
        real dkr = dkx * xr + dky * yr + dkz * zr;
        real pir = pix * xr + piy * yr + piz * zr;
        real pkr = pkx * xr + pky * yr + pkz * zr;

        real r = REAL_SQRT(r2);

        // if use_thole
        real scale3 = uscale[k];
        real scale5 = scale3;
        real damp = pdi * pdamp[k];
        if (damp != 0) {
          real pgamma = REAL_MIN(pti, thole[k]);
          damp = -pgamma * REAL_CUBE(r / damp);
          if (damp > -50) {
            real expdamp = REAL_EXP(damp);
            scale3 *= (1 - expdamp);
            scale5 *= (1 - expdamp * (1 - damp));
          }
        }

        real rr1 = REAL_RECIP(r);
        real rr2 = rr1 * rr1;
        real rr3 = -scale3 * rr1 * rr2;
        real rr5 = 3 * scale5 * rr1 * rr2 * rr2;

        real fid1 = rr3 * dkx + rr5 * dkr * xr;
        real fid2 = rr3 * dky + rr5 * dkr * yr;
        real fid3 = rr3 * dkz + rr5 * dkr * zr;
        real fkd1 = rr3 * dix + rr5 * dir * xr;
        real fkd2 = rr3 * diy + rr5 * dir * yr;
        real fkd3 = rr3 * diz + rr5 * dir * zr;
        real fip1 = rr3 * pkx + rr5 * pkr * xr;
        real fip2 = rr3 * pky + rr5 * pkr * yr;
        real fip3 = rr3 * pkz + rr5 * pkr * zr;
        real fkp1 = rr3 * pix + rr5 * pir * xr;
        real fkp2 = rr3 * piy + rr5 * pir * yr;
        real fkp3 = rr3 * piz + rr5 * pir * zr;

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
TINKER_NAMESPACE_END
