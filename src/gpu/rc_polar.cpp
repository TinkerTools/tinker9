#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_polar.h"
#include "rc_cudart.h"

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

real (*dir_fieldd)[3];
real (*dir_fieldp)[3];

int use_epolar() { return potent::use_polar; }

void get_epolar_type(int& typ, std::string& typ_str) {
  if (limits::use_ewald) {
    typ = elec_ewald;
    typ_str = "EWALD";
  } else {
    typ = elec_ewald;
    typ_str = "COULOMB";
  }
}

void e_polar_data(int op) {
  if (!use_epolar())
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(polarity));
    check_cudart(cudaFree(thole));
    check_cudart(cudaFree(pdamp));
    check_cudart(cudaFree(polarity_inv));

    check_cudart(cudaFree(ep));
    check_cudart(cudaFree(nep));
    check_cudart(cudaFree(vir_ep));

    check_cudart(cudaFree(dir_fieldd));
    check_cudart(cudaFree(dir_fieldp));
  }

  if (op == op_create) {
    get_epolar_type(epolar_electyp, epolar_electyp_str);

    const size_t rs = sizeof(real);
    size_t size;

    check_cudart(cudaMalloc(&polarity, n * rs));
    check_cudart(cudaMalloc(&thole, n * rs));
    check_cudart(cudaMalloc(&pdamp, rs * n));
    check_cudart(cudaMalloc(&polarity_inv, rs * n));
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

    check_cudart(cudaMalloc(&ep, rs));
    check_cudart(cudaMalloc(&nep, sizeof(int)));
    check_cudart(cudaMalloc(&vir_ep, 9 * rs));

    check_cudart(cudaMalloc(&dir_fieldd, 3 * n * rs));
    check_cudart(cudaMalloc(&dir_fieldp, 3 * n * rs));
  }
}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
#define TINKER_GPU_EPOLAR_DEF_(ver)                                            \
  void tinker_gpu_epolar##ver() {                                              \
    if (gpu::epolar_electyp == gpu::elec_coulomb) {                            \
      tinker_gpu_epolar_coulomb##ver();                                        \
    } else if (gpu::epolar_electyp == gpu::elec_ewald) {                       \
      tinker_gpu_epolar_ewald##ver();                                          \
    }                                                                          \
  }
TINKER_GPU_EPOLAR_DEF_(0);
TINKER_GPU_EPOLAR_DEF_(1);
TINKER_GPU_EPOLAR_DEF_(3);
TINKER_GPU_EPOLAR_DEF_(4);
TINKER_GPU_EPOLAR_DEF_(5);
TINKER_GPU_EPOLAR_DEF_(6);
#undef TINKER_GPU_EPOLAR_DEF_
}
