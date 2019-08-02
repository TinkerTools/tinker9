#include "mod_couple.h"
#include "mod_md.h"
#include "util_array.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
void couple_data(rc_op op) {
  if (op & rc_dealloc) {
    check_rt(cudaFree(coupl_obj_.n12));
    check_rt(cudaFree(coupl_obj_.n13));
    check_rt(cudaFree(coupl_obj_.n14));
    check_rt(cudaFree(coupl_obj_.n15));
    check_rt(cudaFree(coupl_obj_.i12));
    check_rt(cudaFree(coupl_obj_.i13));
    check_rt(cudaFree(coupl_obj_.i14));
    check_rt(cudaFree(coupl_obj_.i15));
    check_rt(cudaFree(coupl));
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    check_rt(cudaMalloc(&coupl_obj_.n12, size));
    check_rt(cudaMalloc(&coupl_obj_.n13, size));
    check_rt(cudaMalloc(&coupl_obj_.n14, size));
    check_rt(cudaMalloc(&coupl_obj_.n15, size));
    size = couple_t::maxn12 * n * rs;
    check_rt(cudaMalloc(&coupl_obj_.i12, size));
    size = couple_t::maxn13 * n * rs;
    check_rt(cudaMalloc(&coupl_obj_.i13, size));
    size = couple_t::maxn14 * n * rs;
    check_rt(cudaMalloc(&coupl_obj_.i14, size));
    size = couple_t::maxn15 * n * rs;
    check_rt(cudaMalloc(&coupl_obj_.i15, size));

    size = sizeof(couple_t);
    check_rt(cudaMalloc(&coupl, size));
  }

  if (op & rc_init) {
    size_t size;

    // see also attach.f
    const int maxn13 = 3 * sizes::maxval;
    const int maxn14 = 9 * sizes::maxval;
    const int maxn15 = 27 * sizes::maxval;
    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    size = couple_t::maxn12 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n12[i];
      nbuf[i] = nn;
      int base = i * couple_t::maxn12;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i12[i][j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n12, nbuf.data(), n);
    copyin_array(&coupl_obj_.i12[0][0], ibuf.data(), size);

    size = couple_t::maxn13 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n13[i];
      nbuf[i] = nn;
      int base = i * couple_t::maxn13;
      int bask = i * maxn13;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i13[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n13, nbuf.data(), n);
    copyin_array(&coupl_obj_.i13[0][0], ibuf.data(), size);

    size = couple_t::maxn14 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n14[i];
      nbuf[i] = nn;
      int base = i * couple_t::maxn14;
      int bask = i * maxn14;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i14[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n14, nbuf.data(), n);
    copyin_array(&coupl_obj_.i14[0][0], ibuf.data(), size);

    size = couple_t::maxn15 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n15[i];
      nbuf[i] = nn;
      int base = i * couple_t::maxn15;
      int bask = i * maxn15;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i15[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n15, nbuf.data(), n);
    copyin_array(&coupl_obj_.i15[0][0], ibuf.data(), size);

    size = sizeof(couple_t);
    check_rt(cudaMemcpy(coupl, &coupl_obj_, size, cudaMemcpyHostToDevice));
  }
}
TINKER_NAMESPACE_END
