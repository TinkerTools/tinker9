#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/e_mpole.h"
#include "rc_cudart.h"
#include "util/fort_str.h"

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
real* torque;

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
    check_cudart(cudaFree(torque));
  }

  if (op == op_create) {
    get_empole_type(electyp, electyp_str);

    const size_t rs = sizeof(real);
    size_t size;

    size = sizeof(local_frame_def_st);
    check_cudart(cudaMalloc(&zaxis, n * size));
    std::vector<int> zaxisbuf(4 * n);
    for (int i = 0; i < n; ++i) {
      int base = 4 * i;
      zaxisbuf[base] = mpole::zaxis[i] - 1;
      zaxisbuf[base + 1] = mpole::xaxis[i] - 1;
      zaxisbuf[base + 2] = mpole::yaxis[i] - 1;
      fstr_view str = mpole::polaxe[i];
      int val;
      if (str == "Z-Only")
        val = pole_z_only;
      else if (str == "Z-then-X")
        val = pole_z_then_x;
      else if (str == "Bisector")
        val = pole_bisector;
      else if (str == "Z-Bisect")
        val = pole_z_bisect;
      else if (str == "3-Fold")
        val = pole_3_fold;
      else
        val = pole_none;
      zaxisbuf[base + 3] = val;
    }
    copyin_data(reinterpret_cast<int*>(zaxis), zaxisbuf.data(), 4 * n);

    size = rs * mpl_total;
    check_cudart(cudaMalloc(&pole, n * size));
    check_cudart(cudaMalloc(&rpole, n * size));
    std::vector<double> polebuf(mpl_total * n);
    for (int i = 0; i < n; ++i) {
      int b1 = mpl_total * i;
      int b2 = mpole::maxpole * i;
      for (int j = 0; j < mpl_total; ++j)
        polebuf[b1 + j] = mpole::pole[b2 + mpl_tinker[j]];
    }
    copyin_data(reinterpret_cast<real*>(pole), polebuf.data(), mpl_total * n);

    check_cudart(cudaMalloc(&em, rs));
    check_cudart(cudaMalloc(&nem, sizeof(int)));
    check_cudart(cudaMalloc(&vir_em, 9 * rs));
    if (use_data & use_grad)
      check_cudart(cudaMalloc(&torque, 3 * rs * n));
  }
}
}
TINKER_NAMESPACE_END
