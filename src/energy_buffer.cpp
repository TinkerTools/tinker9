#include "energy_buffer.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
int get_count(Count ne_) {
  int c;
  ne_->sum(&c);
  return c;
}

double get_energy(Energy e_) {
  real real_out;
  e_->sum(&real_out);
  return real_out;
}

void get_virial(double* v_out, Virial vir_) {
  real r6[6];
  vir_->sum(r6);
  for (int i = 0; i < 6; ++i)
    v_out[i] = r6[i];
}

int BondedEnergy::buffer_size() const { return bufsize_; }

Energy BondedEnergy::e() { return e_; }

Virial BondedEnergy::vir() { return vir_; }

void BondedEnergy::dealloc() {
  vir_.close();
  e_.close();
  bufsize_ = 0;
}

void BondedEnergy::alloc(int bsize) {
  if (rc_flag & calc::analyz) {
    e_ = Energy::open();
    vir_ = Virial::open();
  } else {
    e_ = esum_handle;
    vir_ = vir_handle;
  }
  e_->alloc(bsize);
  vir_->alloc(bsize);

  auto esize = e_->size();
  auto vsize = vir_->size();
  assert(esize == vsize);
  bufsize_ = esize;
}

Count NonbondedEnergy::ne() { return ne_; }

void NonbondedEnergy::dealloc() {
  ne_.close();
  BondedEnergy::dealloc();
}

void NonbondedEnergy::alloc(int bsize) {
  BondedEnergy::alloc(bsize);
  ne_ = Count::open();
  ne_->alloc(bsize);

  assert(bufsize_ == ne_->size());
}
TINKER_NAMESPACE_END
