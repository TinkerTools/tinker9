#include "array.h"
#include "md.h"
#include "polgrp.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
static_assert(PolarGroup::maxp11 >= polgrp::maxp11, "");
static_assert(PolarGroup::maxp12 >= polgrp::maxp12, "");
static_assert(PolarGroup::maxp13 >= polgrp::maxp13, "");
static_assert(PolarGroup::maxp14 >= polgrp::maxp14, "");

void polargroup_data(rc_op op) {
  if (op & rc_dealloc) {
    auto& polargroup_obj = polargroup_unit.obj();
    auto* polargroup = polargroup_unit.deviceptr();

    dealloc_bytes(polargroup_obj.np11);
    dealloc_bytes(polargroup_obj.np12);
    dealloc_bytes(polargroup_obj.np13);
    dealloc_bytes(polargroup_obj.np14);
    dealloc_bytes(polargroup_obj.ip11);
    dealloc_bytes(polargroup_obj.ip12);
    dealloc_bytes(polargroup_obj.ip13);
    dealloc_bytes(polargroup_obj.ip14);

    dealloc_bytes(polargroup);

    PolarGroupUnit::clear();
  }

  if (op & rc_alloc) {
    assert(PolarGroupUnit::size() == 0);
    polargroup_unit = PolarGroupUnit::add_new();
    PolarGroup& polargroup_obj = polargroup_unit.obj();
    PolarGroup*& polargroup = polargroup_unit.deviceptr();

    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    alloc_bytes(&polargroup_obj.np11, size);
    alloc_bytes(&polargroup_obj.np12, size);
    alloc_bytes(&polargroup_obj.np13, size);
    alloc_bytes(&polargroup_obj.np14, size);
    size = PolarGroup::maxp11 * n * rs;
    alloc_bytes(&polargroup_obj.ip11, size);
    size = PolarGroup::maxp12 * n * rs;
    alloc_bytes(&polargroup_obj.ip12, size);
    size = PolarGroup::maxp13 * n * rs;
    alloc_bytes(&polargroup_obj.ip13, size);
    size = PolarGroup::maxp14 * n * rs;
    alloc_bytes(&polargroup_obj.ip14, size);

    size = sizeof(PolarGroup);
    alloc_bytes(&polargroup, size);
  }

  if (op & rc_init) {
    auto& polargroup_obj = polargroup_unit.obj();
    auto* polargroup = polargroup_unit.deviceptr();

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
    copyin_array(polargroup_obj.np11, nbuf.data(), n);
    copyin_array(&polargroup_obj.ip11[0][0], ibuf.data(), size);

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
    copyin_array(polargroup_obj.np12, nbuf.data(), n);
    copyin_array(&polargroup_obj.ip12[0][0], ibuf.data(), size);

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
    copyin_array(polargroup_obj.np13, nbuf.data(), n);
    copyin_array(&polargroup_obj.ip13[0][0], ibuf.data(), size);

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
    copyin_array(polargroup_obj.np14, nbuf.data(), n);
    copyin_array(&polargroup_obj.ip14[0][0], ibuf.data(), size);

    size = sizeof(PolarGroup);
    copyin_bytes(polargroup, &polargroup_obj, size);
  }
}
TINKER_NAMESPACE_END
