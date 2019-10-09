#include "couple.h"
#include "dev_array.h"
#include "ext/tinker/detail/couple.hh"
#include "ext/tinker/detail/sizes.hh"
#include "md.h"

TINKER_NAMESPACE_BEGIN
Couple::~Couple() {
  device_array::deallocate(n12, n13, n14, n15, i12, i13, i14, i15);
}

void couple_data(rc_op op) {
  if (op & rc_dealloc)
    CoupleUnit::clear();

  if (op & rc_alloc) {
    assert(CoupleUnit::size() == 0);
    couple_unit = CoupleUnit::open();
    auto& coupl_obj = *couple_unit;

    device_array::allocate(n, //
                           &coupl_obj.n12, &coupl_obj.n13, &coupl_obj.n14,
                           &coupl_obj.n15, //
                           &coupl_obj.i12, &coupl_obj.i13, &coupl_obj.i14,
                           &coupl_obj.i15);
  }

  if (op & rc_init) {
    auto& coupl_obj = *couple_unit;

    // see also attach.f
    const int maxn13 = 3 * sizes::maxval;
    const int maxn14 = 9 * sizes::maxval;
    const int maxn15 = 27 * sizes::maxval;
    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    ibuf.resize(Couple::maxn12 * n);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n12[i];
      nbuf[i] = nn;
      int base = i * Couple::maxn12;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i12[i][j];
        ibuf[base + j] = k - 1;
      }
    }
    device_array::copyin(n, coupl_obj.n12, nbuf.data());
    device_array::copyin(n, coupl_obj.i12, ibuf.data());

    ibuf.resize(Couple::maxn13 * n);
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
    device_array::copyin(n, coupl_obj.n13, nbuf.data());
    device_array::copyin(n, coupl_obj.i13, ibuf.data());

    ibuf.resize(Couple::maxn14 * n);
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
    device_array::copyin(n, coupl_obj.n14, nbuf.data());
    device_array::copyin(n, coupl_obj.i14, ibuf.data());

    ibuf.resize(Couple::maxn15 * n);
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
    device_array::copyin(n, coupl_obj.n15, nbuf.data());
    device_array::copyin(n, coupl_obj.i15, ibuf.data());

    couple_unit.update_deviceptr(coupl_obj);
  }
}
TINKER_NAMESPACE_END
