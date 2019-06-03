#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_polar.h"
#include "rc_cudart.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
real* polarity;
real* pdamp;
real* polarity_inv;

real* ep;
int* nep;
real* vir_ep;

int use_epolar() { return potent::use_polar; }

void e_polar_data(int op) {
  if (!use_epolar())
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(polarity));
    check_cudart(cudaFree(pdamp));
    check_cudart(cudaFree(polarity_inv));

    check_cudart(cudaFree(ep));
    check_cudart(cudaFree(nep));
    check_cudart(cudaFree(vir_ep));
  }

  if (op == op_create) {
    const size_t rs = sizeof(real);
    size_t size;

    // see also polmin in induce.f

    const double polmin = 0.00000001;
    std::vector<double> pbuf(n), pdbuf(n), pinvbuf(n);
    for (int i = 0; i < n; ++i) {
      pbuf[i] = polar::polarity[i];
      pdbuf[i] = polar::pdamp[i];
      pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
    }
    copyin_data(polarity, pbuf.data(), n);
    copyin_data(pdamp, pdbuf.data(), n);
    copyin_data(polarity_inv, pinvbuf.data(), n);
  }
}
}
TINKER_NAMESPACE_END
