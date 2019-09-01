#include "e_tors.h"

#include "ext/tinker/detail/torpot.hh"
#include "ext/tinker/detail/tors.hh"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void etors_data(rc_op op) {
  if (!use_potent(torsion_term))
    return;

  if (op & rc_dealloc) {
    et_handle.dealloc();
  }

  if (op & rc_alloc) {
    ntors = count_bonded_term(torsion_term);
    itors_vec.reserve(4 * ntors);
    tors1_vec.reserve(4 * ntors);
    tors2_vec.reserve(4 * ntors);
    tors3_vec.reserve(4 * ntors);
    tors4_vec.reserve(4 * ntors);
    tors5_vec.reserve(4 * ntors);
    tors6_vec.reserve(4 * ntors);

    et_handle.alloc(ntors);
  }

  if (op & rc_init) {
    std::vector<int> ibuf(4 * ntors);
    for (int i = 0; i < 4 * ntors; ++i) {
      ibuf[i] = tors::itors[i] - 1;
    }
    itors_vec.copyin(ibuf.data(), 4 * ntors);
    tors1_vec.copyin(tors::tors1, 4 * ntors);
    tors2_vec.copyin(tors::tors2, 4 * ntors);
    tors3_vec.copyin(tors::tors3, 4 * ntors);
    tors4_vec.copyin(tors::tors4, 4 * ntors);
    tors5_vec.copyin(tors::tors5, 4 * ntors);
    tors6_vec.copyin(tors::tors6, 4 * ntors);
    torsunit = torpot::torsunit;
  }
}

extern void etors_acc_impl_(int vers);
void etors(int vers) { etors_acc_impl_(vers); }
TINKER_NAMESPACE_END
