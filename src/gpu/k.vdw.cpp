#include "gpu/e.vdw.h"
#include "gpu/mdstate.h"
#include "tinker.mod.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int vdwtyp = 0;
std::string vdwtyp_str;

real *xred, *yred, *zred;

real* ev;
int use_evdw() { return potent::use_vdw; }

void get_evdw_type(int& typ, std::string& typ_str) {
  fstr_view str = vdwpot::vdwtyp;
  typ_str = str.trim();
  if (str == "LENNARD-JONES")
    typ = evdw_lj;
  else if (str == "BUCKINGHAM")
    typ = evdw_buck;
  else if (str == "MM3-HBOND")
    typ = evdw_mm3hb;
  else if (str == "BUFFERED-14-7")
    typ = evdw_hal;
  else if (str == "GAUSSIAN")
    typ = evdw_gauss;
}

real get_evdw() {
  real e;
  check_cudart(cudaMemcpy(&e, ev, sizeof(real), cudaMemcpyDeviceToHost));
  return e;
}

void e_vdw_data(int op) {
  if (!use_evdw())
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(xred));
    check_cudart(cudaFree(yred));
    check_cudart(cudaFree(zred));

    check_cudart(cudaFree(ev));
  }

  if (op == op_create) {
    get_evdw_type(vdwtyp, vdwtyp_str);

    const size_t rs = sizeof(real);
    size_t size;

    size = n * rs;
    check_cudart(cudaMalloc(&xred, size));
    check_cudart(cudaMalloc(&yred, size));
    check_cudart(cudaMalloc(&zred, size));

    check_cudart(cudaMalloc(&ev, rs));
  }
}
}
TINKER_NAMESPACE_END
