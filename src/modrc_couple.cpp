#include "array.h"
#include "couple.h"
#include "md.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
void couple_data(rc_op op) {
  if (op & rc_dealloc) {
    dealloc_bytes(coupl_obj_.n12);
    dealloc_bytes(coupl_obj_.n13);
    dealloc_bytes(coupl_obj_.n14);
    dealloc_bytes(coupl_obj_.n15);
    dealloc_bytes(coupl_obj_.i12);
    dealloc_bytes(coupl_obj_.i13);
    dealloc_bytes(coupl_obj_.i14);
    dealloc_bytes(coupl_obj_.i15);
    dealloc_bytes(coupl);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    alloc_bytes(&coupl_obj_.n12, size);
    alloc_bytes(&coupl_obj_.n13, size);
    alloc_bytes(&coupl_obj_.n14, size);
    alloc_bytes(&coupl_obj_.n15, size);
    size = Couple::maxn12 * n * rs;
    alloc_bytes(&coupl_obj_.i12, size);
    size = Couple::maxn13 * n * rs;
    alloc_bytes(&coupl_obj_.i13, size);
    size = Couple::maxn14 * n * rs;
    alloc_bytes(&coupl_obj_.i14, size);
    size = Couple::maxn15 * n * rs;
    alloc_bytes(&coupl_obj_.i15, size);

    size = sizeof(Couple);
    alloc_bytes(&coupl, size);
  }

  if (op & rc_init) {
    size_t size;

    // see also attach.f
    const int maxn13 = 3 * sizes::maxval;
    const int maxn14 = 9 * sizes::maxval;
    const int maxn15 = 27 * sizes::maxval;
    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    size = Couple::maxn12 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n12[i];
      nbuf[i] = nn;
      int base = i * Couple::maxn12;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i12[i][j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n12, nbuf.data(), n);
    copyin_array(&coupl_obj_.i12[0][0], ibuf.data(), size);

    size = Couple::maxn13 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n13[i];
      nbuf[i] = nn;
      int base = i * Couple::maxn13;
      int bask = i * maxn13;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i13[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n13, nbuf.data(), n);
    copyin_array(&coupl_obj_.i13[0][0], ibuf.data(), size);

    size = Couple::maxn14 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n14[i];
      nbuf[i] = nn;
      int base = i * Couple::maxn14;
      int bask = i * maxn14;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i14[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n14, nbuf.data(), n);
    copyin_array(&coupl_obj_.i14[0][0], ibuf.data(), size);

    size = Couple::maxn15 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n15[i];
      nbuf[i] = nn;
      int base = i * Couple::maxn15;
      int bask = i * maxn15;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i15[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_array(coupl_obj_.n15, nbuf.data(), n);
    copyin_array(&coupl_obj_.i15[0][0], ibuf.data(), size);

    size = sizeof(Couple);
    copyin_bytes(coupl, &coupl_obj_, size);
  }
}
TINKER_NAMESPACE_END
