#include "e_mpole.h"
#include "ext/tinker/detail/couple.hh"
#include "ext/tinker/detail/mplpot.hh"
#include "ext/tinker/detail/sizes.hh"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include "switch.h"

TINKER_NAMESPACE_BEGIN
void empole_data(rc_op op) {
  if (!use_potent(mpole_term))
    return;

  if (op & rc_dealloc) {
    nmpole_excluded_ = 0;

    device_array::deallocate(mpole_excluded_, mpole_excluded_scale_);

    em_handle.dealloc();
  }

  if (op & rc_alloc) {
    m2scale = mplpot::m2scale;
    m3scale = mplpot::m3scale;
    m4scale = mplpot::m4scale;
    m5scale = mplpot::m5scale;

    mpole_excluded_ = 0;
    std::vector<int> exclik;
    std::vector<real> excls;
    // see also attach.f
    const int maxn13 = 3 * sizes::maxval;
    const int maxn14 = 9 * sizes::maxval;
    const int maxn15 = 27 * sizes::maxval;
    for (int i = 0; i < n; ++i) {
      int nn;
      int bask;

      if (m2scale != 1) {
        nn = couple::n12[i];
        for (int j = 0; j < nn; ++j) {
          int k = couple::i12[i][j];
          k -= 1;
          if (k > i) {
            exclik.push_back(i);
            exclik.push_back(k);
            excls.push_back(m2scale - 1);
          }
        }
      }

      if (m3scale != 1) {
        nn = couple::n13[i];
        bask = i * maxn13;
        for (int j = 0; j < nn; ++j) {
          int k = couple::i13[bask + j];
          k -= 1;
          if (k > i) {
            exclik.push_back(i);
            exclik.push_back(k);
            excls.push_back(m3scale - 1);
          }
        }
      }

      if (m4scale != 1) {
        nn = couple::n14[i];
        bask = i * maxn14;
        for (int j = 0; j < nn; ++j) {
          int k = couple::i14[bask + j];
          k -= 1;
          if (k > i) {
            exclik.push_back(i);
            exclik.push_back(k);
            excls.push_back(m4scale - 1);
          }
        }
      }

      if (m5scale != 1) {
        nn = couple::n15[i];
        bask = i * maxn15;
        for (int j = 0; j < nn; ++j) {
          int k = couple::i15[bask + j];
          k -= 1;
          if (k > i) {
            exclik.push_back(i);
            exclik.push_back(k);
            excls.push_back(m5scale - 1);
          }
        }
      }
    }
    nmpole_excluded_ = excls.size();
    device_array::allocate(nmpole_excluded_, &mpole_excluded_,
                           &mpole_excluded_scale_);
    device_array::copyin(nmpole_excluded_, mpole_excluded_, exclik.data());
    device_array::copyin(nmpole_excluded_, mpole_excluded_scale_, excls.data());

    em_handle.alloc(n);
  }

  if (op & rc_init) {
    if (use_ewald()) {
      empole_electyp = elec_t::ewald;
    } else {
      empole_electyp = elec_t::coulomb;
    }

    if (empole_electyp == elec_t::coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);
  }
}

void empole(int vers) {
  if (empole_electyp == elec_t::coulomb)
    empole_coulomb(vers);
  else if (empole_electyp == elec_t::ewald)
    empole_ewald(vers);
}
TINKER_NAMESPACE_END
