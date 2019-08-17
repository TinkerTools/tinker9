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
  real r9[9];
  vir_->sum(r9);
  for (int i = 0; i < 9; ++i)
    v_out[i] = r9[i];
}

int BondedEnergy::buffer_size() const { return bufsize_; }

Energy BondedEnergy::e() { return e_; }

Virial BondedEnergy::vir() { return vir_; }

void BondedEnergy::dealloc() { bufsize_ = 0; }

void BondedEnergy::alloc(int bsize) {
  if (rc_flag & calc::analyz) {
    e_ = Energy::inquire();
    vir_ = Virial::inquire();
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

void NonbondedEnergy::dealloc() { bufsize_ = 0; }

void NonbondedEnergy::alloc(int bsize) {
  BondedEnergy::alloc(bsize);
  ne_ = Count::inquire();
  ne_->alloc(bsize);

  assert(bufsize_ == ne_->size());
}
TINKER_NAMESPACE_END
