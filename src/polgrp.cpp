#include "polgrp.h"
#include "dev_array.h"
#include "ext/tinker/detail/polgrp.hh"
#include "md.h"

TINKER_NAMESPACE_BEGIN
static_assert(PolarGroup::maxp11 >= polgrp::maxp11, "");
static_assert(PolarGroup::maxp12 >= polgrp::maxp12, "");
static_assert(PolarGroup::maxp13 >= polgrp::maxp13, "");
static_assert(PolarGroup::maxp14 >= polgrp::maxp14, "");

PolarGroup::~PolarGroup() {
  device_array::deallocate(np11, np12, np13, np14, ip11, ip12, ip13, ip14);
}

void polargroup_data(rc_op op) {
  if (op & rc_dealloc)
    PolarGroupUnit::clear();

  if (op & rc_alloc) {
    assert(PolarGroupUnit::size() == 0);
    polargroup_unit = PolarGroupUnit::open();
    auto& polargroup_obj = *polargroup_unit;

    device_array::allocate(n, //
                           &polargroup_obj.np11, &polargroup_obj.np12,
                           &polargroup_obj.np13,
                           &polargroup_obj.np14, //
                           &polargroup_obj.ip11, &polargroup_obj.ip12,
                           &polargroup_obj.ip13, &polargroup_obj.ip14);
  }

  if (op & rc_init) {
    auto& polargroup_obj = *polargroup_unit;

    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    ibuf.resize(PolarGroup::maxp11 * n);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np11[i];
      nbuf[i] = nn;
      int base = i * PolarGroup::maxp11;
      int bask = i * polgrp::maxp11;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip11[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    device_array::copyin(n, polargroup_obj.np11, nbuf.data());
    device_array::copyin(n, polargroup_obj.ip11, ibuf.data());

    ibuf.resize(PolarGroup::maxp12 * n);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np12[i];
      nbuf[i] = nn;
      int base = i * PolarGroup::maxp12;
      int bask = i * polgrp::maxp12;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip12[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    device_array::copyin(n, polargroup_obj.np12, nbuf.data());
    device_array::copyin(n, polargroup_obj.ip12, ibuf.data());

    ibuf.resize(PolarGroup::maxp13 * n);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np13[i];
      nbuf[i] = nn;
      int base = i * PolarGroup::maxp13;
      int bask = i * polgrp::maxp13;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip13[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    device_array::copyin(n, polargroup_obj.np13, nbuf.data());
    device_array::copyin(n, polargroup_obj.ip13, ibuf.data());

    ibuf.resize(PolarGroup::maxp14 * n);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np14[i];
      nbuf[i] = nn;
      int base = i * PolarGroup::maxp14;
      int bask = i * polgrp::maxp14;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip14[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    device_array::copyin(n, polargroup_obj.np14, nbuf.data());
    device_array::copyin(n, polargroup_obj.ip14, ibuf.data());

    polargroup_unit.init_deviceptr(polargroup_obj);
  }
}
TINKER_NAMESPACE_END
