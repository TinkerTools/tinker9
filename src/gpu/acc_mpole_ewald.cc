#include "acc_e.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void empole_real_self_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  const real f = chgpot::electric / chgpot::dielec;

  const real off = ewald_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_obj_.maxnlst;

  const real m2scale = mplpot::m2scale;
  const real m3scale = mplpot::m3scale;
  const real m4scale = mplpot::m4scale;
  const real m5scale = mplpot::m5scale;

  static std::vector<real> mscalebuf;
  mscalebuf.resize(n, 1);
  real* mscale = mscalebuf.data();

  const int pu = epme_unit;
  const real aewald = pme_obj(pu).aewald;
  const real aewald_sq_2 = 2 * aewald * aewald;
  const real fterm = -f * aewald * 0.5 * M_2_SQRTPI;

  real bn[6];

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,box,couple,mlst,\
                        rpole,\
                        em,nem,vir_em,trqx,trqy,trqz)\
              firstprivate(mscale[0:n])
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int n12i = couple->n12[i];
    const int n13i = couple->n13[i];
    const int n14i = couple->n14[i];
    const int n15i = couple->n15[i];
    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      mscale[couple->i12[i][j]] = m2scale;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      mscale[couple->i13[i][j]] = m3scale;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      mscale[couple->i14[i][j]] = m4scale;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      mscale[couple->i15[i][j]] = m5scale;

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

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent private(bn[0:6])
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        real r = REAL_SQRT(r2);
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
        real dik = dix * dkx + diy * dky + diz * dkz;
        real qik = qix * qkx + qiy * qky + qiz * qkz;
        real diqk = dix * qkx + diy * qky + diz * qkz;
        real dkqi = dkx * qix + dky * qiy + dkz * qiz;
        real qiqk = 2 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) +
            qixx * qkxx + qiyy * qkyy + qizz * qkzz;

        real invr1 = REAL_RECIP(r);
        real rr2 = invr1 * invr1;

        real rr1 = f * invr1;
        real rr3 = rr1 * rr2;
        real rr5 = 3 * rr3 * rr2;
        real rr7 = 5 * rr5 * rr2;
        real rr9 = 7 * rr7 * rr2;
        MAYBE_UNUSED real rr11;

        real ralpha = aewald * r;

        bn[0] = REAL_ERFC(ralpha) * invr1;
        real alsq2 = 2 * REAL_SQ(aewald);
        real alsq2n = (aewald > 0 ? REAL_RECIP(sqrtpi * aewald) : 0);
        real exp2a = REAL_EXP(-REAL_SQ(ralpha));
        if_constexpr(!do_g) {
          #pragma acc loop seq
          for (int j = 1; j <= 4; ++j) {
            alsq2n *= alsq2;
            bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
          }
        }
        else {
          #pragma acc loop seq
          for (int j = 1; j <= 5; ++j) {
            alsq2n *= alsq2;
            bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
          }
        }
        bn[0] *= f;
        bn[1] *= f;
        bn[2] *= f;
        bn[3] *= f;
        bn[4] *= f;
        if_constexpr(do_g) {
          bn[5] *= f;
          rr11 = 9 * rr9 * rr2;
        }

        real term1 = ci * ck;
        real term2 = ck * dir - ci * dkr + dik;
        real term3 = ci * qkr + ck * qir - dir * dkr + 2 * (dkqi - diqk + qiqk);
        real term4 = dir * qkr - dkr * qir - 4 * qik;
        real term5 = qir * qkr;

        if_constexpr(do_a) {
          real efull = term1 * rr1 + term2 * rr3 + term3 * rr5 + term4 * rr7 +
              term5 * rr9;
          efull *= mscale[k];
          if (efull != 0) {
            #pragma acc atomic update
            *nem += 1;
          }
        } // end if (do_a)

        real scalek = 1 - mscale[k];
        rr1 = bn[0] - scalek * rr1;
        rr3 = bn[1] - scalek * rr3;
        rr5 = bn[2] - scalek * rr5;
        rr7 = bn[3] - scalek * rr7;
        rr9 = bn[4] - scalek * rr9;
        if_constexpr(do_e) {
          real e = term1 * rr1 + term2 * rr3 + term3 * rr5 + term4 * rr7 +
              term5 * rr9;
          #pragma acc atomic update
          *em += e;
        } // end if (do_e)

        if_constexpr(do_g) {

          // gradient

          real qixk = qixx * qkx + qixy * qky + qixz * qkz;
          real qiyk = qixy * qkx + qiyy * qky + qiyz * qkz;
          real qizk = qixz * qkx + qiyz * qky + qizz * qkz;
          real qkxi = qkxx * qix + qkxy * qiy + qkxz * qiz;
          real qkyi = qkxy * qix + qkyy * qiy + qkyz * qiz;
          real qkzi = qkxz * qix + qkyz * qiy + qkzz * qiz;

          real diqkx = dix * qkxx + diy * qkxy + diz * qkxz;
          real diqky = dix * qkxy + diy * qkyy + diz * qkyz;
          real diqkz = dix * qkxz + diy * qkyz + diz * qkzz;
          real dkqix = dkx * qixx + dky * qixy + dkz * qixz;
          real dkqiy = dkx * qixy + dky * qiyy + dkz * qiyz;
          real dkqiz = dkx * qixz + dky * qiyz + dkz * qizz;

          rr11 = bn[5] - scalek * rr11;

          real de = term1 * rr3 + term2 * rr5 + term3 * rr7 + term4 * rr9 +
              term5 * rr11;

          term1 = -ck * rr3 + dkr * rr5 - qkr * rr7;
          term2 = ci * rr3 + dir * rr5 + qir * rr7;
          term3 = 2 * rr5;
          term4 = 2 * (-ck * rr5 + dkr * rr7 - qkr * rr9);
          term5 = 2 * (-ci * rr5 - dir * rr7 - qir * rr9);
          real term6 = 4 * rr7;

          real frcx = de * xr + term1 * dix + term2 * dkx +
              term3 * (diqkx - dkqix) + term4 * qix + term5 * qkx +
              term6 * (qixk + qkxi);
          real frcy = de * yr + term1 * diy + term2 * dky +
              term3 * (diqky - dkqiy) + term4 * qiy + term5 * qky +
              term6 * (qiyk + qkyi);
          real frcz = de * zr + term1 * diz + term2 * dkz +
              term3 * (diqkz - dkqiz) + term4 * qiz + term5 * qkz +
              term6 * (qizk + qkzi);

          #pragma acc atomic update
          gx[i] += frcx;
          #pragma acc atomic update
          gy[i] += frcy;
          #pragma acc atomic update
          gz[i] += frcz;
          #pragma acc atomic update
          gx[k] -= frcx;
          #pragma acc atomic update
          gy[k] -= frcy;
          #pragma acc atomic update
          gz[k] -= frcz;

          // torque

          real dirx = diy * zr - diz * yr;
          real diry = diz * xr - dix * zr;
          real dirz = dix * yr - diy * xr;
          real dkrx = dky * zr - dkz * yr;
          real dkry = dkz * xr - dkx * zr;
          real dkrz = dkx * yr - dky * xr;
          real dikx = diy * dkz - diz * dky;
          real diky = diz * dkx - dix * dkz;
          real dikz = dix * dky - diy * dkx;

          real qirx = qiz * yr - qiy * zr;
          real qiry = qix * zr - qiz * xr;
          real qirz = qiy * xr - qix * yr;
          real qkrx = qkz * yr - qky * zr;
          real qkry = qkx * zr - qkz * xr;
          real qkrz = qky * xr - qkx * yr;
          real qikx = qky * qiz - qkz * qiy;
          real qiky = qkz * qix - qkx * qiz;
          real qikz = qkx * qiy - qky * qix;

          real qikrx = qizk * yr - qiyk * zr;
          real qikry = qixk * zr - qizk * xr;
          real qikrz = qiyk * xr - qixk * yr;
          real qkirx = qkzi * yr - qkyi * zr;
          real qkiry = qkxi * zr - qkzi * xr;
          real qkirz = qkyi * xr - qkxi * yr;

          real diqkrx = diqkz * yr - diqky * zr;
          real diqkry = diqkx * zr - diqkz * xr;
          real diqkrz = diqky * xr - diqkx * yr;
          real dkqirx = dkqiz * yr - dkqiy * zr;
          real dkqiry = dkqix * zr - dkqiz * xr;
          real dkqirz = dkqiy * xr - dkqix * yr;

          real dqikx = diy * qkz - diz * qky + dky * qiz - dkz * qiy -
              2 *
                  (qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy -
                   qiyz * qkyy - qizz * qkyz);
          real dqiky = diz * qkx - dix * qkz + dkz * qix - dkx * qiz -
              2 *
                  (qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz -
                   qixy * qkyz - qixz * qkzz);
          real dqikz = dix * qky - diy * qkx + dkx * qiy - dky * qix -
              2 *
                  (qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx -
                   qiyy * qkxy - qiyz * qkxz);

          real ttmi[3], ttmk[3];
          ttmi[0] = -rr3 * dikx + term1 * dirx + term3 * (dqikx + dkqirx) -
              term4 * qirx - term6 * (qikrx + qikx);
          ttmi[1] = -rr3 * diky + term1 * diry + term3 * (dqiky + dkqiry) -
              term4 * qiry - term6 * (qikry + qiky);
          ttmi[2] = -rr3 * dikz + term1 * dirz + term3 * (dqikz + dkqirz) -
              term4 * qirz - term6 * (qikrz + qikz);
          ttmk[0] = rr3 * dikx + term2 * dkrx - term3 * (dqikx + diqkrx) -
              term5 * qkrx - term6 * (qkirx - qikx);
          ttmk[1] = rr3 * diky + term2 * dkry - term3 * (dqiky + diqkry) -
              term5 * qkry - term6 * (qkiry - qiky);
          ttmk[2] = rr3 * dikz + term2 * dkrz - term3 * (dqikz + diqkrz) -
              term5 * qkrz - term6 * (qkirz - qikz);

          #pragma acc atomic update
          trqx[i] += ttmi[0];
          #pragma acc atomic update
          trqy[i] += ttmi[1];
          #pragma acc atomic update
          trqz[i] += ttmi[2];
          #pragma acc atomic update
          trqx[k] += ttmk[0];
          #pragma acc atomic update
          trqy[k] += ttmk[1];
          #pragma acc atomic update
          trqz[k] += ttmk[2];

          // virial

          if_constexpr(do_v) {
            real vxx = -xr * frcx;
            real vxy = -0.5f * (yr * frcx + xr * frcy);
            real vxz = -0.5f * (zr * frcx + xr * frcz);
            real vyy = -yr * frcy;
            real vyz = -0.5f * (zr * frcy + yr * frcz);
            real vzz = -zr * frcz;

            #pragma acc atomic update
            vir_em[_xx] += vxx;
            #pragma acc atomic update
            vir_em[_yx] += vxy;
            #pragma acc atomic update
            vir_em[_zx] += vxz;
            #pragma acc atomic update
            vir_em[_xy] += vxy;
            #pragma acc atomic update
            vir_em[_yy] += vyy;
            #pragma acc atomic update
            vir_em[_zy] += vyz;
            #pragma acc atomic update
            vir_em[_xz] += vxz;
            #pragma acc atomic update
            vir_em[_yz] += vyz;
            #pragma acc atomic update
            vir_em[_zz] += vzz;
          } // end if (do_v)
        }   // end if (do_g)
      }     // end if (r2 <= off2)
    }       // end for (int kk)

    // reset exclusion coefficients for connected atoms

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      mscale[couple->i12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      mscale[couple->i13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      mscale[couple->i14[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      mscale[couple->i15[i][j]] = 1;

    // compute the self-energy part of the Ewald summation

    real cii = ci * ci;
    real dii = dix * dix + diy * diy + diz * diz;
    real qii = 2 * (qixy * qixy + qixz * qixz + qiyz * qiyz) + qixx * qixx +
        qiyy * qiyy + qizz * qizz;

    if_constexpr(do_e) {
      real e = fterm *
          (cii + aewald_sq_2 * (dii / 3 + 2 * aewald_sq_2 * qii * (real)0.2));
      #pragma acc atomic update
      *em += e;
      if_constexpr(do_a) {
        #pragma acc atomic update
        *nem += 1;
      }
    } // end if (do_e)
  }   // end for (int i)
}

template <int USE>
void empole_recip_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  const int pu = epme_unit;
  cmp_to_fmp(pu, cmp, fmp);
  grid_mpole(pu, fmp);
  fftfront(pu);
  if_constexpr(do_v) {
    if (vir_m) {
      pme_conv1(pu, vir_m);
      #pragma acc parallel loop independent deviceptr(vir_m,vir_em)
      for (int i = 0; i < 9; ++i) {
        vir_em[i] += vir_m[i];
      }
    } else {
      pme_conv1(pu, vir_em);
    }
  }
  else {
    pme_conv0(pu);
  }
  fftback(pu);
  fphi_mpole(pu, fphi);
  fphi_to_cphi(pu, fphi, cphi);

  constexpr int deriv1[] = {2, 5, 8, 9, 11, 16, 18, 14, 15, 20};
  constexpr int deriv2[] = {3, 8, 6, 10, 14, 12, 19, 16, 20, 17};
  constexpr int deriv3[] = {4, 9, 10, 7, 15, 17, 13, 20, 18, 19};

  pme_st& st = pme_obj(pu);
  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const real f = chgpot::electric / chgpot::dielec;

  #pragma acc parallel loop independent\
              deviceptr(gx,gy,gz,box,\
                        cmp,fmp,cphi,fphi,\
                        em,vir_em,trqx,trqy,trqz)
  for (int i = 0; i < n; ++i) {
    real e = 0;
    real f1 = 0;
    real f2 = 0;
    real f3 = 0;
    #pragma acc loop independent reduction(+:e,f1,f2,f3)
    for (int k = 0; k < 10; ++k) {
      if_constexpr(do_e) { e += fmp[i][k] * fphi[i][k]; }

      if_constexpr(do_g) {
        f1 += fmp[i][k] * fphi[i][deriv1[k] - 1];
        f2 += fmp[i][k] * fphi[i][deriv2[k] - 1];
        f3 += fmp[i][k] * fphi[i][deriv3[k] - 1];
      }
    } // end for (int k)

    // increment the permanent multipole energy and gradient

    if_constexpr(do_e) {
      #pragma acc atomic update
      *em += 0.5f * e * f;
    }

    if_constexpr(do_g) {
      f1 *= nfft1;
      f2 *= nfft2;
      f3 *= nfft3;

      real h1 =
          box->recip[0][0] * f1 + box->recip[1][0] * f2 + box->recip[2][0] * f3;
      real h2 =
          box->recip[0][1] * f1 + box->recip[1][1] * f2 + box->recip[2][1] * f3;
      real h3 =
          box->recip[0][2] * f1 + box->recip[1][2] * f2 + box->recip[2][2] * f3;

      #pragma acc atomic update
      gx[i] += h1 * f;
      #pragma acc atomic update
      gy[i] += h2 * f;
      #pragma acc atomic update
      gz[i] += h3 * f;

      // resolve site torques then increment forces and virial

      real tem1 = cmp[i][3] * cphi[i][2] - cmp[i][2] * cphi[i][3] +
          2 * (cmp[i][6] - cmp[i][5]) * cphi[i][9] + cmp[i][8] * cphi[i][7] +
          cmp[i][9] * cphi[i][5] - cmp[i][7] * cphi[i][8] -
          cmp[i][9] * cphi[i][6];
      real tem2 = cmp[i][1] * cphi[i][3] - cmp[i][3] * cphi[i][1] +
          2 * (cmp[i][4] - cmp[i][6]) * cphi[i][8] + cmp[i][7] * cphi[i][9] +
          cmp[i][8] * cphi[i][6] - cmp[i][8] * cphi[i][4] -
          cmp[i][9] * cphi[i][7];
      real tem3 = cmp[i][2] * cphi[i][1] - cmp[i][1] * cphi[i][2] +
          2 * (cmp[i][5] - cmp[i][4]) * cphi[i][7] + cmp[i][7] * cphi[i][4] +
          cmp[i][9] * cphi[i][8] - cmp[i][7] * cphi[i][5] -
          cmp[i][8] * cphi[i][9];
      tem1 *= f;
      tem2 *= f;
      tem3 *= f;

      #pragma acc atomic update
      trqx[i] += tem1;
      #pragma acc atomic update
      trqy[i] += tem2;
      #pragma acc atomic update
      trqz[i] += tem3;

      if_constexpr(do_v) {
        real vxx = -cmp[i][1] * cphi[i][1] - 2 * cmp[i][4] * cphi[i][4] -
            cmp[i][7] * cphi[i][7] - cmp[i][8] * cphi[i][8];
        real vxy = -0.5f * (cmp[i][2] * cphi[i][1] + cmp[i][1] * cphi[i][2]) -
            (cmp[i][4] + cmp[i][5]) * cphi[i][7] -
            0.5f * cmp[i][7] * (cphi[i][4] + cphi[i][5]) -
            0.5f * (cmp[i][8] * cphi[i][9] + cmp[i][9] * cphi[i][8]);
        real vxz = -0.5f * (cmp[i][3] * cphi[i][1] + cmp[i][1] * cphi[i][3]) -
            (cmp[i][4] + cmp[i][6]) * cphi[i][8] -
            0.5f * cmp[i][8] * (cphi[i][4] + cphi[i][6]) -
            0.5f * (cmp[i][7] * cphi[i][9] + cmp[i][9] * cphi[i][7]);
        real vyy = -cmp[i][2] * cphi[i][2] - 2 * cmp[i][5] * cphi[i][5] -
            cmp[i][7] * cphi[i][7] - cmp[i][9] * cphi[i][9];
        real vyz = -0.5f * (cmp[i][3] * cphi[i][2] + cmp[i][2] * cphi[i][3]) -
            (cmp[i][5] + cmp[i][6]) * cphi[i][9] -
            0.5f * cmp[i][9] * (cphi[i][5] + cphi[i][6]) -
            0.5f * (cmp[i][7] * cphi[i][8] + cmp[i][8] * cphi[i][7]);
        real vzz = -cmp[i][3] * cphi[i][3] - 2 * cmp[i][6] * cphi[i][6] -
            cmp[i][8] * cphi[i][8] - cmp[i][9] * cphi[i][9];
        vxx *= f;
        vxy *= f;
        vxz *= f;
        vyy *= f;
        vyz *= f;
        vzz *= f;

        #pragma acc atomic update
        vir_em[_xx] += vxx;
        #pragma acc atomic update
        vir_em[_yx] += vxy;
        #pragma acc atomic update
        vir_em[_zx] += vxz;
        #pragma acc atomic update
        vir_em[_xy] += vxy;
        #pragma acc atomic update
        vir_em[_yy] += vyy;
        #pragma acc atomic update
        vir_em[_zy] += vyz;
        #pragma acc atomic update
        vir_em[_xz] += vxz;
        #pragma acc atomic update
        vir_em[_yz] += vyz;
        #pragma acc atomic update
        vir_em[_zz] += vzz;
      } // end if (do_v)
    }   // end if (do_g)
  }     // end for (int i)
}

template <int USE>
void empole_ewald_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  if_constexpr(do_e || do_a || do_v) {
    #pragma acc serial deviceptr(em,nem,vir_em)
    {
      if_constexpr(do_e) { *em = 0; }
      if_constexpr(do_a) { *nem = 0; }
      if_constexpr(do_v) {
        for (int i = 0; i < 9; ++i)
          vir_em[i] = 0;
      }
    }
  }

  empole_real_self_tmpl<USE>();

  empole_recip_tmpl<USE>();
}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
void tinker_gpu_empole_ewald0() { gpu::empole_ewald_tmpl<gpu::v0>(); }
void tinker_gpu_empole_ewald1() { gpu::empole_ewald_tmpl<gpu::v1>(); }
void tinker_gpu_empole_ewald3() { gpu::empole_ewald_tmpl<gpu::v3>(); }
void tinker_gpu_empole_ewald4() { gpu::empole_ewald_tmpl<gpu::v4>(); }
void tinker_gpu_empole_ewald5() { gpu::empole_ewald_tmpl<gpu::v5>(); }
void tinker_gpu_empole_ewald6() { gpu::empole_ewald_tmpl<gpu::v6>(); }
}
