#include "gpu/decl_mdstate.h"
#include "gpu/rc.h"
#include <map>

TINKER_NAMESPACE_BEGIN
namespace gpu {
real* esum;
real epot, eksum, ekin[3][3];
real *gx, *gy, *gz;
real* vir;

//======================================================================

namespace energy_buffer_ {
static const int virlen = 16;
static int end;
static int cap;
static int* nebuf;
static real* ebuf;
static real* vbuf;
static std::map<int**, int> ne_addr_idx;
static std::map<real**, int> e_addr_idx;
static std::map<real**, int> v_addr_idx;

static void data(rc_t rc) {
  if (rc & rc_dealloc) {
    end = 0;
    cap = 0;

    check_cudart(cudaFree(nebuf));
    nebuf = nullptr;
    check_cudart(cudaFree(ebuf));
    ebuf = nullptr;
    check_cudart(cudaFree(vbuf));
    vbuf = nullptr;

    ne_addr_idx.clear();
    e_addr_idx.clear();
    v_addr_idx.clear();
  }

  if (rc & rc_alloc) {
    end = 0;
    cap = 4; // default initial capacity

    const size_t rs = sizeof(real);
    check_cudart(cudaMalloc(&nebuf, sizeof(int) * cap));
    check_cudart(cudaMalloc(&ebuf, rs * cap));
    check_cudart(cudaMalloc(&vbuf, rs * cap * virlen));
  }

  if (rc & rc_copyin) {
    zero_array(nebuf, cap);
    zero_array(ebuf, cap);
    zero_array(vbuf, cap * virlen);
  }
}

static void grow_if_must() {
  if (end < cap)
    return;

  int old_cap;
  old_cap = cap;
  cap *= 2;

  const size_t rs = sizeof(real);

  int* new_nebuf;
  check_cudart(cudaMalloc(&new_nebuf, sizeof(int) * cap));
  check_cudart(cudaMemcpy(new_nebuf, nebuf, sizeof(int) * old_cap,
                          cudaMemcpyDeviceToDevice));
  check_cudart(cudaFree(nebuf));
  nebuf = new_nebuf;

  real* new_ebuf;
  check_cudart(cudaMalloc(&new_ebuf, rs * cap));
  check_cudart(
      cudaMemcpy(new_ebuf, ebuf, rs * old_cap, cudaMemcpyDeviceToDevice));
  check_cudart(cudaFree(ebuf));
  ebuf = new_ebuf;

  real* new_vbuf;
  check_cudart(cudaMalloc(&new_vbuf, rs * cap * virlen));
  check_cudart(cudaMemcpy(new_vbuf, vbuf, rs * old_cap * virlen,
                          cudaMemcpyDeviceToDevice));
  check_cudart(cudaFree(vbuf));
  vbuf = new_vbuf;

  for (auto it : ne_addr_idx) {
    int** p = it.first;
    int idx = it.second;
    *p = &nebuf[idx];
  }

  for (auto it : e_addr_idx) {
    real** p = it.first;
    int idx = it.second;
    *p = &ebuf[idx];
  }

  for (auto it : v_addr_idx) {
    real** p = it.first;
    int idx = it.second;
    *p = &vbuf[idx * virlen];
  }
}

static void alloc3(real** pe, real** pv, int** pne) {
  auto it = e_addr_idx.find(pe);
  if (it != e_addr_idx.end())
    return;

  grow_if_must();
  *pe = ebuf + end;
  e_addr_idx[pe] = end;
  *pv = vbuf + end * virlen;
  v_addr_idx[pv] = end;
  if (pne) {
    *pne = nebuf + end;
    ne_addr_idx[pne] = end;
  }

  ++end;
}

static void dealloc3(real* /* pe */, real* /* pv */, int* /* pne */) {}
}
//======================================================================

void egv_data(rc_t rc) {
  if (use_data & (use_analyz | use_energy | use_virial)) {
    energy_buffer_::data(rc);

    if (rc & rc_dealloc) {
      free_ev(esum, vir);
    }

    if (rc & rc_alloc) {
      alloc_ev(&esum, &vir);
    }

    if (rc & rc_copyin) {
      if (use_energy & use_data)
        copyin_array(esum, &energi::esum, 1);

      if (use_virial & use_data)
        copyin_array(vir, &virial::vir[0][0], 9);
    }
  }

  if (use_data & use_grad) {
    if (rc & rc_dealloc) {
      check_cudart(cudaFree(gx));
      check_cudart(cudaFree(gy));
      check_cudart(cudaFree(gz));
    }

    if (rc & rc_alloc) {
      const size_t size = sizeof(real) * n;
      check_cudart(cudaMalloc(&gx, size));
      check_cudart(cudaMalloc(&gy, size));
      check_cudart(cudaMalloc(&gz, size));
    }

    // We can never assume whether or not deriv::desum was allocated, because it
    // was allocated inside subroutine gradient(...), which would be skipped in
    // subroutine mdinit() if a dyn file existed to restart a simulation.

    // if (rc & rc_copyin) {
    //   copyin_array2(2, 3, gz, deriv::desum, n);
    //   copyin_array2(2, 3, gz, deriv::desum, n);
    //   copyin_array2(2, 3, gz, deriv::desum, n);
    // }
  }
}

//======================================================================

void alloc_ev(real** gpu_e, real** gpu_v) {
  energy_buffer_::alloc3(gpu_e, gpu_v, nullptr);
}

void free_ev(real* gpu_e, real* gpu_v) {
  energy_buffer_::dealloc3(gpu_e, gpu_v, nullptr);
}

void alloc_nev(int** gpu_ne, real** gpu_e, real** gpu_v) {
  energy_buffer_::alloc3(gpu_e, gpu_v, gpu_ne);
}

void free_nev(int* gpu_ne, real* gpu_e, real* gpu_v) {
  energy_buffer_::dealloc3(gpu_e, gpu_v, gpu_ne);
}

//======================================================================

double get_energy(const real* e_gpu) {
  double e_out;
  copyout_array(&e_out, e_gpu, 1);
  return e_out;
}

int get_count(const int* ecount_gpu) {
  int c;
  copyout_array(&c, ecount_gpu, 1);
  return c;
}

void get_virial(double* v_out, const real* v_gpu) {
  copyout_array(v_out, v_gpu, 9);
}

void zero_egv(int vers) {
  if (vers & use_analyz)
    zero_array(energy_buffer_::nebuf, energy_buffer_::cap);

  if (vers & use_energy)
    zero_array(energy_buffer_::ebuf, energy_buffer_::cap);

  if (vers & use_virial)
    zero_array(energy_buffer_::vbuf,
               energy_buffer_::cap * energy_buffer_::virlen);

  if (vers & use_grad) {
    zero_array(gx, n);
    zero_array(gy, n);
    zero_array(gz, n);
  }
}

void zero_egv() { zero_egv(use_data & vmask); }

extern void sum_energy_acc_impl_(real* ebuf, int end);
extern void sum_virial_acc_impl_(real* vbuf, int end, int virlen);
void sum_energies(int vers) {
  if (vers & use_energy)
    sum_energy_acc_impl_(energy_buffer_::ebuf, energy_buffer_::end);

  if (vers & use_virial)
    sum_virial_acc_impl_(energy_buffer_::vbuf, energy_buffer_::end,
                         energy_buffer_::virlen);
}
}
TINKER_NAMESPACE_END
