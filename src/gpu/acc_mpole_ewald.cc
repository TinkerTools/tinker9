#include "acc_e.h"
#include "gpu/acc_fmat.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"
#include <vector>

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void empole_real_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  const real f = chgpot::electric / chgpot::dielec;
  const real aewald = pme_obj(epme_unit).aewald;
  const real aewald_sq_2 = 2 * aewald * aewald;
  const real fterm = -f * aewald * 0.5 * M_2_SQRTPI;

  #pragma acc data deviceptr(x,y,z,gx,gy,gz,box,couple,mlst,\
                             rpole,\
                             em,nem,vir_em,trqx,trqy,trqz)
  {
    #pragma acc parallel loop
    for (int i = 0; i < n; ++i) {

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
      } // end if_constexpr(do_e)
    }
  }
}

template <int USE>
void empole_recip_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  cmp_to_fmp(fmp, epme_unit);
  grid_mpole(fmp);
  fftfront(epme_unit);
  if_constexpr(do_v) { pme_conv1(epme_unit, vir_em); }
  else {
    pme_conv0(epme_unit);
  }
  fftback(epme_unit);
  fphi_mpole(fphi);
  fphi_to_cphi(fphi, cphi, epme_unit);

  constexpr int deriv1[] = {2, 5, 8, 9, 11, 16, 18, 14, 15, 20};
  constexpr int deriv2[] = {3, 8, 6, 10, 14, 12, 19, 16, 20, 17};
  constexpr int deriv3[] = {4, 9, 10, 7, 15, 17, 13, 20, 18, 19};

  pme_st& st = pme_obj(epme_unit);
  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const real f = chgpot::electric / chgpot::dielec;

  #pragma acc kernels deviceptr(gx,gy,gz,box,\
                                rpole,fmp,cphi,fphi,\
                                em,vir_em,trqx,trqy,trqz)
  {
    #pragma acc loop independent
    for (int iatom = 0; iatom < n; ++iatom) {
      real e = 0;
      real f1 = 0;
      real f2 = 0;
      real f3 = 0;
      #pragma acc loop independent reduction(+:e,f1,f2,f3)
      for (int k = 0; k < 10; ++k) {
        if_constexpr(do_e) { e += fmp[iatom][k] * fphi[iatom][k]; }

        if_constexpr(do_g) {
          f1 += fmp[iatom][k] * fphi[iatom][deriv1[k] - 1];
          f2 += fmp[iatom][k] * fphi[iatom][deriv2[k] - 1];
          f3 += fmp[iatom][k] * fphi[iatom][deriv3[k] - 1];
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

        fmat_real3 recip(box->recip);
        real h1 = recip(1, 1) * f1 + recip(1, 2) * f2 + recip(1, 3) * f3;
        real h2 = recip(2, 1) * f1 + recip(2, 2) * f2 + recip(2, 3) * f3;
        real h3 = recip(3, 1) * f1 + recip(3, 2) * f2 + recip(3, 3) * f3;

        #pragma acc atomic update
        gx[iatom] += h1 * f;
        #pragma acc atomic update
        gy[iatom] += h2 * f;
        #pragma acc atomic update
        gz[iatom] += h3 * f;

        // resolve site torques then increment forces and virial

        real cmp_cpp[10];
        farray_real cmp(cmp_cpp);
        fmat_real<10> cphi(gpu::cphi);

        const int i = iatom + 1;

        cmp(1) = rpole[iatom][mpl_pme_0];
        cmp(2) = rpole[iatom][mpl_pme_x];
        cmp(3) = rpole[iatom][mpl_pme_y];
        cmp(4) = rpole[iatom][mpl_pme_z];
        cmp(5) = rpole[iatom][mpl_pme_xx];
        cmp(6) = rpole[iatom][mpl_pme_yy];
        cmp(7) = rpole[iatom][mpl_pme_zz];
        cmp(8) = 2 * rpole[iatom][mpl_pme_xy];
        cmp(9) = 2 * rpole[iatom][mpl_pme_xz];
        cmp(10) = 2 * rpole[iatom][mpl_pme_yz];

        real tem1 = cmp(4) * cphi(3, i) - cmp(3) * cphi(4, i) +
            2 * (cmp(7) - cmp(6)) * cphi(10, i) + cmp(9) * cphi(8, i) +
            cmp(10) * cphi(6, i) - cmp(8) * cphi(9, i) - cmp(10) * cphi(7, i);
        real tem2 = cmp(2) * cphi(4, i) - cmp(4) * cphi(2, i) +
            2 * (cmp(5) - cmp(7)) * cphi(9, i) + cmp(8) * cphi(10, i) +
            cmp(9) * cphi(7, i) - cmp(9) * cphi(5, i) - cmp(10) * cphi(8, i);
        real tem3 = cmp(3) * cphi(2, i) - cmp(2) * cphi(3, i) +
            2 * (cmp(6) - cmp(5)) * cphi(8, i) + cmp(8) * cphi(5, i) +
            cmp(10) * cphi(9, i) - cmp(8) * cphi(6, i) - cmp(9) * cphi(10, i);

        #pragma acc atomic update
        trqx[iatom] += tem1;
        #pragma acc atomic update
        trqy[iatom] += tem2;
        #pragma acc atomic update
        trqz[iatom] += tem3;

        if_constexpr(do_v) {
          real vxx = -cmp(2) * cphi(2, i) - 2 * cmp(5) * cphi(5, i) -
              cmp(8) * cphi(8, i) - cmp(9) * cphi(9, i);
          real vxy = -0.5f * (cmp(3) * cphi(2, i) + cmp(2) * cphi(3, i)) -
              (cmp(5) + cmp(6)) * cphi(8, i) -
              0.5f * cmp(8) * (cphi(5, i) + cphi(6, i)) -
              0.5f * (cmp(9) * cphi(10, i) + cmp(10) * cphi(9, i));
          real vxz = -0.5f * (cmp(4) * cphi(2, i) + cmp(2) * cphi(4, i)) -
              (cmp(5) + cmp(7)) * cphi(9, i) -
              0.5f * cmp(9) * (cphi(5, i) + cphi(7, i)) -
              0.5f * (cmp(8) * cphi(10, i) + cmp(10) * cphi(8, i));
          real vyy = -cmp(3) * cphi(3, i) - 2 * cmp(6) * cphi(6, i) -
              cmp(8) * cphi(8, i) - cmp(10) * cphi(10, i);
          real vyz = -0.5f * (cmp(4) * cphi(3, i) + cmp(3) * cphi(4, i)) -
              (cmp(6) + cmp(7)) * cphi(10, i) -
              0.5f * cmp(10) * (cphi(6, i) + cphi(7, i)) -
              0.5f * (cmp(8) * cphi(9, i) + cmp(9) * cphi(8, i));
          real vzz = -cmp(4) * cphi(4, i) - 2 * cmp(7) * cphi(7, i) -
              cmp(9) * cphi(9, i) - cmp(10) * cphi(10, i);

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
}

template <int USE>
void empole_ewald_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  zero_torque();

  chkpole();
  rotpole();

  #pragma acc data deviceptr(em,nem,vir_em)
  #pragma acc serial
  {
    *em = 0;
    if_constexpr(do_a) { *nem = 0; }
    if_constexpr(do_v) {
      for (int i = 0; i < 9; ++i)
        vir_em[i] = 0;
    }
  }

  empole_real_tmpl<USE>();

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
