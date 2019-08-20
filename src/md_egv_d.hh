#include "array.h"
#include "ext/tinker/detail/energi.hh"
#include "ext/tinker/detail/virial.hh"
#include "md.h"
#include <map>

TINKER_NAMESPACE_BEGIN
namespace {
static bool use_ev_() {
  return rc_flag & (calc::analyz | calc::energy | calc::virial);
}

template <class T>
struct Buffer {
  static const int virlen = 16;
  int end;
  int cap;
  int* nebuf;
  T* ebuf;
  T* vbuf;
  std::map<int**, int> ne_addr_idx;
  std::map<T**, int> e_addr_idx;
  std::map<T**, int> v_addr_idx;

  void data(rc_op op) {
    if (!use_ev_())
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

      const size_t rs = sizeof(T);
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

  void grow_if_must() {
    if (end < cap)
      return;

    int old_cap;
    old_cap = cap;
    cap *= 2;

    const size_t rs = sizeof(T);

    int* new_nebuf;
    alloc_bytes(&new_nebuf, sizeof(int) * cap);
    copy_bytes(new_nebuf, nebuf, sizeof(int) * old_cap);
    dealloc_bytes(nebuf);
    nebuf = new_nebuf;

    T* new_ebuf;
    alloc_bytes(&new_ebuf, rs * cap);
    copy_bytes(new_ebuf, ebuf, rs * old_cap);
    dealloc_bytes(ebuf);
    ebuf = new_ebuf;

    T* new_vbuf;
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
      T** p = it.first;
      int idx = it.second;
      *p = &ebuf[idx];
    }

    for (auto it : v_addr_idx) {
      T** p = it.first;
      int idx = it.second;
      *p = &vbuf[idx * virlen];
    }
  }

  void alloc3(T** pe, T** pv, int** pne) {
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

  void dealloc3(T*,   // pe
                T*,   // pv
                int*) // pne
  {}
};

static Buffer<double> double_buffer;
static void double_data(rc_op op) { return double_buffer.data(op); }
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
static void ev_data_(rc_op op) {
  if (!use_ev_())
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
    const size_t size = sizeof(double) * n;
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
  rc_man energy_double42_{double_data, op};
  rc_man ev42_{ev_data_, op};
  rc_man grad42_{grad_data_, op};
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void alloc_ev(real_buffer_t** gpu_e, real_buffer_t** gpu_v) {
  double_buffer.alloc3(gpu_e, gpu_v, nullptr);
}

void dealloc_ev(real_buffer_t* gpu_e, real_buffer_t* gpu_v) {
  double_buffer.dealloc3(gpu_e, gpu_v, nullptr);
}

void alloc_nev(int** gpu_ne, real_buffer_t** gpu_e, real_buffer_t** gpu_v) {
  double_buffer.alloc3(gpu_e, gpu_v, gpu_ne);
}

void dealloc_nev(int* gpu_ne, real_buffer_t* gpu_e, real_buffer_t* gpu_v) {
  double_buffer.dealloc3(gpu_e, gpu_v, gpu_ne);
}

double get_energy(const real_buffer_t* e_gpu) {
  double e_out;
  copyout_array(&e_out, e_gpu, 1);
  return e_out;
}

int get_count(const int* ecount_gpu) {
  int c;
  copyout_array(&c, ecount_gpu, 1);
  return c;
}

void get_virial(double* v_out, const real_buffer_t* v_gpu) {
  copyout_array(v_out, v_gpu, 9);
}

void zero_egv(int vers) {
  if (vers & calc::analyz) {
    zero_array(double_buffer.nebuf, double_buffer.cap);
  }

  if (vers & calc::energy) {
    zero_array(double_buffer.ebuf, double_buffer.cap);
  }

  if (vers & calc::virial) {
    zero_array(double_buffer.vbuf, double_buffer.cap * double_buffer.virlen);
  }

  if (vers & calc::grad) {
    zero_array(gx, n);
    zero_array(gy, n);
    zero_array(gz, n);
  }
}

void zero_egv() { zero_egv(rc_flag & calc::vmask); }

extern void sum_energy_acc_impl_(double* ebuf, int end);
extern void sum_virial_acc_impl_(double* vbuf, int end, int virlen);
void sum_energies(int vers) {
  if (vers & calc::energy)
    sum_energy_acc_impl_(double_buffer.ebuf, double_buffer.end);

  if (vers & calc::virial)
    sum_virial_acc_impl_(double_buffer.vbuf, double_buffer.end,
                         double_buffer.virlen);
}
TINKER_NAMESPACE_END
