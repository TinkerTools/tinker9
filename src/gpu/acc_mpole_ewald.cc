#include "acc_e.h"
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
  // fftfront(epme_unit);
  // if_constexpr(do_v) {
  //   pme_conv1(epme_unit, vir);
  // } else {
  //   pme_conv0(epme_unit);
  // }
  // fftback(epme_unit);
  // fphi_mpole(fphi);
  // fphi_to_cphi(fphi, cphi, epme_unit);
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
