#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/e_mpole.h"
#include "rc_cudart.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int electyp;
std::string electyp_str;

local_frame_def_st* zaxis;

real (*pole)[mpl_total];
real (*rpole)[mpl_total];

real* em;
int* nem;
real* vir_em;
int use_empole() { return potent::use_mpole; }

void get_empole_type(int& typ, std::string& typ_str) {
  if (limits::use_ewald) {
    typ = elec_ewald;
    typ_str = "EWALD";
  } else {
    typ = elec_coulomb;
    typ_str = "COULOMB";
  }
}
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_mlist_build() {
  m_tinker_using_namespace;
  if (gpu::use_mpole_list()) {
    gpu::nblist_construct(gpu::mlist_obj_, gpu::mlst);
  }
}
void tinker_gpu_mlist_update() {
  m_tinker_using_namespace;
  if (gpu::use_mpole_list()) {
    gpu::nblist_update(gpu::mlist_obj_, gpu::mlst);
  }
}
}

TINKER_NAMESPACE_BEGIN
namespace gpu {
void e_mpole_data(int op) {
  if (!use_empole())
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(zaxis));
    check_cudart(cudaFree(pole));
    check_cudart(cudaFree(rpole));

    check_cudart(cudaFree(em));
    check_cudart(cudaFree(nem));
    check_cudart(cudaFree(vir_em));
  }

  if (op == op_create) {
    get_empole_type(electyp, electyp_str);

    const size_t rs = sizeof(real);
    size_t size;

    size = sizeof(local_frame_def_st);
    check_cudart(cudaMalloc(&zaxis, n * size));
    size = rs * mpl_total;
    check_cudart(cudaMalloc(&pole, n * size));
    check_cudart(cudaMalloc(&rpole, n * size));

    check_cudart(cudaMalloc(&em, rs));
    check_cudart(cudaMalloc(&nem, sizeof(int)));
    check_cudart(cudaMalloc(&vir_em, 9 * rs));
  }
}
}
TINKER_NAMESPACE_END
