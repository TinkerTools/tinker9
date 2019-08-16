#include "couple.h"
#include "array.h"
#include "ext/tinker/detail/couple.hh"
#include "ext/tinker/detail/sizes.hh"
#include "md.h"

TINKER_NAMESPACE_BEGIN
Couple::~Couple() {
  dealloc_bytes(n12);
  dealloc_bytes(n13);
  dealloc_bytes(n14);
  dealloc_bytes(n15);
  dealloc_bytes(i12);
  dealloc_bytes(i13);
  dealloc_bytes(i14);
  dealloc_bytes(i15);
}

void couple_data(rc_op op) {
  if (op & rc_dealloc)
    CoupleUnit::clear();

  if (op & rc_alloc) {
    assert(CoupleUnit::size() == 0);
    couple_unit = CoupleUnit::alloc_new();
    auto& coupl_obj = *couple_unit;

    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    alloc_bytes(&coupl_obj.n12, size);
    alloc_bytes(&coupl_obj.n13, size);
    alloc_bytes(&coupl_obj.n14, size);
    alloc_bytes(&coupl_obj.n15, size);
    size = Couple::maxn12 * n * rs;
    alloc_bytes(&coupl_obj.i12, size);
    size = Couple::maxn13 * n * rs;
    alloc_bytes(&coupl_obj.i13, size);
    size = Couple::maxn14 * n * rs;
    alloc_bytes(&coupl_obj.i14, size);
    size = Couple::maxn15 * n * rs;
    alloc_bytes(&coupl_obj.i15, size);
  }

  if (op & rc_init) {
    auto& coupl_obj = *couple_unit;
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
    copyin_array(coupl_obj.n12, nbuf.data(), n);
    copyin_array(&coupl_obj.i12[0][0], ibuf.data(), size);

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
    copyin_array(coupl_obj.n13, nbuf.data(), n);
    copyin_array(&coupl_obj.i13[0][0], ibuf.data(), size);

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
    copyin_array(coupl_obj.n14, nbuf.data(), n);
    copyin_array(&coupl_obj.i14[0][0], ibuf.data(), size);

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
    copyin_array(coupl_obj.n15, nbuf.data(), n);
    copyin_array(&coupl_obj.i15[0][0], ibuf.data(), size);

    couple_unit.init_deviceptr(coupl_obj);
  }
}
TINKER_NAMESPACE_END
