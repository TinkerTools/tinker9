#include "polgrp.h"
#include "dev_mem.h"
#include "ext/tinker/detail/polgrp.hh"
#include "md.h"

TINKER_NAMESPACE_BEGIN
static_assert(PolarGroup::maxp11 >= polgrp::maxp11, "");
static_assert(PolarGroup::maxp12 >= polgrp::maxp12, "");
static_assert(PolarGroup::maxp13 >= polgrp::maxp13, "");
static_assert(PolarGroup::maxp14 >= polgrp::maxp14, "");

PolarGroup::~PolarGroup() {
  DeviceMemory::deallocate_bytes(np11);
  DeviceMemory::deallocate_bytes(np12);
  DeviceMemory::deallocate_bytes(np13);
  DeviceMemory::deallocate_bytes(np14);
  DeviceMemory::deallocate_bytes(ip11);
  DeviceMemory::deallocate_bytes(ip12);
  DeviceMemory::deallocate_bytes(ip13);
  DeviceMemory::deallocate_bytes(ip14);
}

void polargroup_data(rc_op op) {
  if (op & rc_dealloc)
    PolarGroupUnit::clear();

  if (op & rc_alloc) {
    assert(PolarGroupUnit::size() == 0);
    polargroup_unit = PolarGroupUnit::open();
    auto& polargroup_obj = *polargroup_unit;
    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    DeviceMemory::allocate_bytes(&polargroup_obj.np11, size);
    DeviceMemory::allocate_bytes(&polargroup_obj.np12, size);
    DeviceMemory::allocate_bytes(&polargroup_obj.np13, size);
    DeviceMemory::allocate_bytes(&polargroup_obj.np14, size);
    size = PolarGroup::maxp11 * n * rs;
    DeviceMemory::allocate_bytes(&polargroup_obj.ip11, size);
    size = PolarGroup::maxp12 * n * rs;
    DeviceMemory::allocate_bytes(&polargroup_obj.ip12, size);
    size = PolarGroup::maxp13 * n * rs;
    DeviceMemory::allocate_bytes(&polargroup_obj.ip13, size);
    size = PolarGroup::maxp14 * n * rs;
    DeviceMemory::allocate_bytes(&polargroup_obj.ip14, size);
  }

  if (op & rc_init) {
    auto& polargroup_obj = *polargroup_unit;
    size_t size;

    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    size = PolarGroup::maxp11 * n;
    ibuf.resize(size);
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
    DeviceMemory::copyin_array(polargroup_obj.np11, nbuf.data(), n);
    DeviceMemory::copyin_array(&polargroup_obj.ip11[0][0], ibuf.data(), size);

    size = PolarGroup::maxp12 * n;
    ibuf.resize(size);
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
    DeviceMemory::copyin_array(polargroup_obj.np12, nbuf.data(), n);
    DeviceMemory::copyin_array(&polargroup_obj.ip12[0][0], ibuf.data(), size);

    size = PolarGroup::maxp13 * n;
    ibuf.resize(size);
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
    DeviceMemory::copyin_array(polargroup_obj.np13, nbuf.data(), n);
    DeviceMemory::copyin_array(&polargroup_obj.ip13[0][0], ibuf.data(), size);

    size = PolarGroup::maxp14 * n;
    ibuf.resize(size);
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
    DeviceMemory::copyin_array(polargroup_obj.np14, nbuf.data(), n);
    DeviceMemory::copyin_array(&polargroup_obj.ip14[0][0], ibuf.data(), size);

    polargroup_unit.init_deviceptr(polargroup_obj);
  }
}
TINKER_NAMESPACE_END
