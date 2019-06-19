#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_polgrp.h"
#include "rc_cudart.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
static_assert(polargroup_st::maxp11 >= polgrp::maxp11, "");
static_assert(polargroup_st::maxp12 >= polgrp::maxp12, "");
static_assert(polargroup_st::maxp13 >= polgrp::maxp13, "");
static_assert(polargroup_st::maxp14 >= polgrp::maxp14, "");

polargroup_st polargroup_obj_;
polargroup_st* polargroup;

void polargroup_data(int op) {
  if (op & op_dealloc) {
    check_cudart(cudaFree(polargroup_obj_.np11));
    check_cudart(cudaFree(polargroup_obj_.np12));
    check_cudart(cudaFree(polargroup_obj_.np13));
    check_cudart(cudaFree(polargroup_obj_.np14));
    check_cudart(cudaFree(polargroup_obj_.ip11));
    check_cudart(cudaFree(polargroup_obj_.ip12));
    check_cudart(cudaFree(polargroup_obj_.ip13));
    check_cudart(cudaFree(polargroup_obj_.ip14));
    check_cudart(cudaFree(polargroup));
  }

  // TODO
  if (op & op_alloc) {
  }
  if (op & op_copyin) {
  }

  if (op == op_create) {
    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    check_cudart(cudaMalloc(&polargroup_obj_.np11, size));
    check_cudart(cudaMalloc(&polargroup_obj_.np12, size));
    check_cudart(cudaMalloc(&polargroup_obj_.np13, size));
    check_cudart(cudaMalloc(&polargroup_obj_.np14, size));
    size = polargroup_st::maxp11 * n * rs;
    check_cudart(cudaMalloc(&polargroup_obj_.ip11, size));
    size = polargroup_st::maxp12 * n * rs;
    check_cudart(cudaMalloc(&polargroup_obj_.ip12, size));
    size = polargroup_st::maxp13 * n * rs;
    check_cudart(cudaMalloc(&polargroup_obj_.ip13, size));
    size = polargroup_st::maxp14 * n * rs;
    check_cudart(cudaMalloc(&polargroup_obj_.ip14, size));

    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    size = polargroup_st::maxp11 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np11[i];
      nbuf[i] = nn;
      int base = i * polargroup_st::maxp11;
      int bask = i * polgrp::maxp11;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip11[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data(polargroup_obj_.np11, nbuf.data(), n);
    copyin_data(&polargroup_obj_.ip11[0][0], ibuf.data(), size);

    size = polargroup_st::maxp12 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np12[i];
      nbuf[i] = nn;
      int base = i * polargroup_st::maxp12;
      int bask = i * polgrp::maxp12;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip12[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data(polargroup_obj_.np12, nbuf.data(), n);
    copyin_data(&polargroup_obj_.ip12[0][0], ibuf.data(), size);

    size = polargroup_st::maxp13 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np13[i];
      nbuf[i] = nn;
      int base = i * polargroup_st::maxp13;
      int bask = i * polgrp::maxp13;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip13[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data(polargroup_obj_.np13, nbuf.data(), n);
    copyin_data(&polargroup_obj_.ip13[0][0], ibuf.data(), size);

    size = polargroup_st::maxp14 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = polgrp::np14[i];
      nbuf[i] = nn;
      int base = i * polargroup_st::maxp14;
      int bask = i * polgrp::maxp14;
      for (int j = 0; j < nn; ++j) {
        int k = polgrp::ip14[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data(polargroup_obj_.np14, nbuf.data(), n);
    copyin_data(&polargroup_obj_.ip14[0][0], ibuf.data(), size);

    size = sizeof(polargroup_st);
    check_cudart(cudaMalloc(&polargroup, size));
    check_cudart(
        cudaMemcpy(polargroup, &polargroup_obj_, size, cudaMemcpyHostToDevice));
  }
}
}
TINKER_NAMESPACE_END
