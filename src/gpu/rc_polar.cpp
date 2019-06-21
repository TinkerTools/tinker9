#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/decl_switch.h"
#include "gpu/e_polar.h"
#include "rc_cudart.h"
#include "util/format_print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int epolar_electyp;
std::string epolar_electyp_str;

real* polarity;
real* thole;
real* pdamp;
real* polarity_inv;

real* ep;
int* nep;
real* vir_ep;

real (*ufld)[3];
real (*dufld)[6];

real (*work01__)[3];
real (*work02__)[3];
real (*work03__)[3];
real (*work04__)[3];
real (*work05__)[3];
real (*work06__)[3];
real (*work07__)[3];
real (*work08__)[3];
real (*work09__)[3];
real (*work10__)[3];

int use_epolar() { return potent::use_polar; }

void get_epolar_type(int& typ, std::string& typ_str) {
  if (use_ewald()) {
    typ = elec_ewald;
    typ_str = "EWALD";
  } else {
    typ = elec_coulomb;
    typ_str = "COULOMB";
  }
}

void epolar_data(int op) {
  if (!use_epolar())
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(polarity));
    check_cudart(cudaFree(thole));
    check_cudart(cudaFree(pdamp));
    check_cudart(cudaFree(polarity_inv));

    check_cudart(cudaFree(ep));
    check_cudart(cudaFree(nep));
    check_cudart(cudaFree(vir_ep));

    check_cudart(cudaFree(ufld));
    check_cudart(cudaFree(dufld));

    check_cudart(cudaFree(work01__));
    check_cudart(cudaFree(work02__));
    check_cudart(cudaFree(work03__));
    check_cudart(cudaFree(work04__));
    check_cudart(cudaFree(work05__));
    check_cudart(cudaFree(work06__));
    check_cudart(cudaFree(work07__));
    check_cudart(cudaFree(work08__));
    check_cudart(cudaFree(work09__));
    check_cudart(cudaFree(work10__));
  }

  if (op & op_alloc) {
    const size_t rs = sizeof(real);
    size_t size;

    check_cudart(cudaMalloc(&polarity, n * rs));
    check_cudart(cudaMalloc(&thole, n * rs));
    check_cudart(cudaMalloc(&pdamp, rs * n));
    check_cudart(cudaMalloc(&polarity_inv, rs * n));

    check_cudart(cudaMalloc(&ep, rs));
    check_cudart(cudaMalloc(&nep, sizeof(int)));
    check_cudart(cudaMalloc(&vir_ep, 9 * rs));

    if (use_data & use_grad) {
      check_cudart(cudaMalloc(&ufld, rs * 3 * n));
      check_cudart(cudaMalloc(&dufld, rs * 6 * n));
    } else {
      ufld = nullptr;
      dufld = nullptr;
    }

    check_cudart(cudaMalloc(&work01__, 3 * n * rs));
    check_cudart(cudaMalloc(&work02__, 3 * n * rs));
    check_cudart(cudaMalloc(&work03__, 3 * n * rs));
    check_cudart(cudaMalloc(&work04__, 3 * n * rs));
    check_cudart(cudaMalloc(&work05__, 3 * n * rs));
    check_cudart(cudaMalloc(&work06__, 3 * n * rs));
    check_cudart(cudaMalloc(&work07__, 3 * n * rs));
    check_cudart(cudaMalloc(&work08__, 3 * n * rs));
    check_cudart(cudaMalloc(&work09__, 3 * n * rs));
    check_cudart(cudaMalloc(&work10__, 3 * n * rs));
  }

  if (op & op_copyin) {
    get_epolar_type(epolar_electyp, epolar_electyp_str);

    if (epolar_electyp == elec_coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);

    // see also polmin in induce.f
    const double polmin = 0.00000001;
    std::vector<double> pinvbuf(n);
    for (int i = 0; i < n; ++i) {
      pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
    }
    copyin_data(polarity, polar::polarity, n);
    copyin_data(thole, polar::thole, n);
    copyin_data(pdamp, polar::pdamp, n);
    copyin_data(polarity_inv, pinvbuf.data(), n);
  }
}

void dfield(real* gpu_field, real* gpu_fieldp) {
  if (epolar_electyp == elec_ewald)
    dfield_ewald(gpu_field, gpu_fieldp);
  else
    dfield_coulomb(gpu_field, gpu_fieldp);
}

void ufield(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
            real* gpu_fieldp) {
  if (epolar_electyp == elec_ewald)
    ufield_ewald(gpu_uind, gpu_uinp, gpu_field, gpu_fieldp);
  else
    ufield_coulomb(gpu_uind, gpu_uinp, gpu_field, gpu_fieldp);
}

void induce(real* gpu_ud, real* gpu_up) {
  induce_mutual_pcg1(gpu_ud, gpu_up);

  if (inform::debug && use_epolar()) {
    std::vector<double> uindbuf;
    uindbuf.resize(3 * n);
    copyout_data(uindbuf.data(), gpu_ud, 3 * n);
    bool header = true;
    for (int i = 0; i < n; ++i) {
      if (polar::polarity[i] != 0) {
        if (header) {
          header = false;
          print(stdout, "\n Induced Dipole Moments (Debye) :\n");
          print(stdout, "\n{0:4s}Atom{0:15s}X{0:12s}Y{0:12s}Z{0:11s}Total\n\n",
                "");
        }
        double u1 = uindbuf[3 * i];
        double u2 = uindbuf[3 * i + 1];
        double u3 = uindbuf[3 * i + 2];
        double unorm = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);
        u1 *= units::debye;
        u2 *= units::debye;
        u3 *= units::debye;
        unorm *= units::debye;
        print(stdout, "{:>8d}     {:13.4f}{:13.4f}{:13.4f} {:13.4f}\n", i + 1,
              u1, u2, u3, unorm);
      }
    }
  }
}

void epolar(int vers) {
  if (epolar_electyp == elec_coulomb)
    epolar_coulomb(vers);
  else if (epolar_electyp == elec_ewald)
    epolar_ewald(vers);
}
}
TINKER_NAMESPACE_END
