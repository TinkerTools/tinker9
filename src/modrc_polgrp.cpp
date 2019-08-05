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
    dealloc_bytes(polargroup_obj_.np11);
    dealloc_bytes(polargroup_obj_.np12);
    dealloc_bytes(polargroup_obj_.np13);
    dealloc_bytes(polargroup_obj_.np14);
    dealloc_bytes(polargroup_obj_.ip11);
    dealloc_bytes(polargroup_obj_.ip12);
    dealloc_bytes(polargroup_obj_.ip13);
    dealloc_bytes(polargroup_obj_.ip14);

    dealloc_bytes(polargroup);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    alloc_bytes(&polargroup_obj_.np11, size);
    alloc_bytes(&polargroup_obj_.np12, size);
    alloc_bytes(&polargroup_obj_.np13, size);
    alloc_bytes(&polargroup_obj_.np14, size);
    size = PolarGroup::maxp11 * n * rs;
    alloc_bytes(&polargroup_obj_.ip11, size);
    size = PolarGroup::maxp12 * n * rs;
    alloc_bytes(&polargroup_obj_.ip12, size);
    size = PolarGroup::maxp13 * n * rs;
    alloc_bytes(&polargroup_obj_.ip13, size);
    size = PolarGroup::maxp14 * n * rs;
    alloc_bytes(&polargroup_obj_.ip14, size);

    size = sizeof(PolarGroup);
    alloc_bytes(&polargroup, size);
  }

  if (op & rc_init) {
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
    copyin_array(polargroup_obj_.np11, nbuf.data(), n);
    copyin_array(&polargroup_obj_.ip11[0][0], ibuf.data(), size);

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
    copyin_array(polargroup_obj_.np12, nbuf.data(), n);
    copyin_array(&polargroup_obj_.ip12[0][0], ibuf.data(), size);

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
    copyin_array(polargroup_obj_.np13, nbuf.data(), n);
    copyin_array(&polargroup_obj_.ip13[0][0], ibuf.data(), size);

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
    copyin_array(polargroup_obj_.np14, nbuf.data(), n);
    copyin_array(&polargroup_obj_.ip14[0][0], ibuf.data(), size);

    size = sizeof(PolarGroup);
    copyin_bytes(polargroup, &polargroup_obj_, size);
  }
}
TINKER_NAMESPACE_END
