#include "array.h"
#include "ext/tinker/detail/energi.hh"
#include "ext/tinker/detail/virial.hh"
#include "md.h"
#include <map>

TINKER_NAMESPACE_BEGIN
namespace energy_buffer_ {
static bool use_ev() {
  return rc_flag & (calc::analyz | calc::energy | calc::virial);
}

static const int virlen = 16;
static int end;
static int cap;
static int* nebuf;
static real* ebuf;
static real* vbuf;
static std::map<int**, int> ne_addr_idx;
static std::map<real**, int> e_addr_idx;
static std::map<real**, int> v_addr_idx;

static void data(rc_op op) {
  if (!use_ev())
    return;

  if (op & rc_dealloc) {
    end = 0;
    cap = 0;

    dealloc_bytes(nebuf);
    nebuf = nullptr;
    dealloc_bytes(ebuf);
    ebuf = nullptr;
    dealloc_bytes(vbuf);
    vbuf = nullptr;

    ne_addr_idx.clear();
    e_addr_idx.clear();
    v_addr_idx.clear();
  }

  if (op & rc_alloc) {
    end = 0;
    cap = 4; // default initial capacity

    const size_t rs = sizeof(real);
    alloc_bytes(&nebuf, sizeof(int) * cap);
    alloc_bytes(&ebuf, rs * cap);
    alloc_bytes(&vbuf, rs * cap * virlen);
  }

  if (op & rc_init) {
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
  alloc_bytes(&new_nebuf, sizeof(int) * cap);
  copy_bytes(new_nebuf, nebuf, sizeof(int) * old_cap);
  dealloc_bytes(nebuf);
  nebuf = new_nebuf;

  real* new_ebuf;
  alloc_bytes(&new_ebuf, rs * cap);
  copy_bytes(new_ebuf, ebuf, rs * old_cap);
  dealloc_bytes(ebuf);
  ebuf = new_ebuf;

  real* new_vbuf;
  alloc_bytes(&new_vbuf, rs * cap * virlen);
  copy_bytes(new_vbuf, vbuf, rs * old_cap * virlen);
  dealloc_bytes(vbuf);
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

static void dealloc3(real*, // pe
                     real*, // pv
                     int*)  // pne
{}
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
static void ev_data_(rc_op op) {
  if (!energy_buffer_::use_ev())
    return;

  if (op & rc_dealloc)
    dealloc_ev(esum, vir);

  if (op & rc_alloc)
    alloc_ev(&esum, &vir);

  if (op & rc_init) {
    if (calc::energy & rc_flag)
      copyin_array(esum, &energi::esum, 1);

    if (calc::virial & rc_flag)
      copyin_array(vir, &virial::vir[0][0], 9);
  }
}

static void grad_data_(rc_op op) {
  if (!(rc_flag & calc::grad))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(gx);
    dealloc_bytes(gy);
    dealloc_bytes(gz);
  }

  if (op & rc_alloc) {
    const size_t size = sizeof(real) * n;
    alloc_bytes(&gx, size);
    alloc_bytes(&gy, size);
    alloc_bytes(&gz, size);
  }

  // we can never assume whether or not deriv::desum was allocated, because it
  // was allocated inside subroutine gradient(...), which would be skipped in
  // subroutine mdinit() if a dyn file existed to restart a simulation.

  // if (op & rc_init) {
  // copy in deriv::sum to gx, gy, and gz
  // }
}

void egv_data(rc_op op) {
  rc_man energy_buffer42_{energy_buffer_::data, op};
  rc_man ev42_{ev_data_, op};
  rc_man grad42_{grad_data_, op};
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void alloc_ev(real** gpu_e, real** gpu_v) {
  energy_buffer_::alloc3(gpu_e, gpu_v, nullptr);
}

void dealloc_ev(real* gpu_e, real* gpu_v) {
  energy_buffer_::dealloc3(gpu_e, gpu_v, nullptr);
}

void alloc_nev(int** gpu_ne, real** gpu_e, real** gpu_v) {
  energy_buffer_::alloc3(gpu_e, gpu_v, gpu_ne);
}

void dealloc_nev(int* gpu_ne, real* gpu_e, real* gpu_v) {
  energy_buffer_::dealloc3(gpu_e, gpu_v, gpu_ne);
}

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
  if (vers & calc::analyz)
    zero_array(energy_buffer_::nebuf, energy_buffer_::cap);

  if (vers & calc::energy)
    zero_array(energy_buffer_::ebuf, energy_buffer_::cap);

  if (vers & calc::virial)
    zero_array(energy_buffer_::vbuf,
               energy_buffer_::cap * energy_buffer_::virlen);

  if (vers & calc::grad) {
    zero_array(gx, n);
    zero_array(gy, n);
    zero_array(gz, n);
  }
}

void zero_egv() { zero_egv(rc_flag & calc::vmask); }

extern void sum_energy_acc_impl_(real* ebuf, int end);
extern void sum_virial_acc_impl_(real* vbuf, int end, int virlen);
void sum_energies(int vers) {
  if (vers & calc::energy)
    sum_energy_acc_impl_(energy_buffer_::ebuf, energy_buffer_::end);

  if (vers & calc::virial)
    sum_virial_acc_impl_(energy_buffer_::vbuf, energy_buffer_::end,
                         energy_buffer_::virlen);
}
TINKER_NAMESPACE_END
