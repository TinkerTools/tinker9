#include "e_tortor.h"
#include "array.h"
#include "ext/tinker/detail/atomid.hh"
#include "ext/tinker/detail/atoms.hh"
#include "ext/tinker/detail/bitor.hh"
#include "ext/tinker/detail/couple.hh"
#include "ext/tinker/detail/torpot.hh"
#include "ext/tinker/detail/tortor.hh"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void etortor_data(rc_op op) {
  if (!use_potent(tortor_term))
    return;

  if (op & rc_dealloc) {
    ett_handle.dealloc();
  }

  if (op & rc_alloc) {
    nbitor = bitor_::nbitor;
    ibitor_vec.reserve(5 * nbitor);

    itt_vec.reserve(3 * nbitor);

    tnx_vec.reserve(ktrtor::maxntt);
    tny_vec.reserve(ktrtor::maxntt);
    int count = ktrtor::maxtgrd * ktrtor::maxntt;
    ttx_vec.reserve(count);
    tty_vec.reserve(count);
    count = ktrtor::maxtgrd2 * ktrtor::maxntt;
    tbf_vec.reserve(count);
    tbx_vec.reserve(count);
    tby_vec.reserve(count);
    tbxy_vec.reserve(count);

    ntortor = count_bonded_term(tortor_term);
    chkttor_ia_vec_.reserve(ntortor);

    ett_handle.alloc(ntortor);
  }

  if (op & rc_init) {
    std::vector<int> ibuf;

    ibuf.resize(5 * nbitor);
    for (int i = 0; i < 5 * nbitor; ++i)
      ibuf[i] = bitor_::ibitor[i] - 1;
    ibitor_vec.copyin(ibuf.data(), 5 * nbitor);

    ibuf.resize(3 * nbitor);
    for (int i = 0; i < 3 * nbitor; ++i)
      ibuf[i] = tortor::itt[i] - 1;
    itt_vec.copyin(ibuf.data(), 3 * nbitor);

    ibuf.resize(ktrtor::maxntt);
    for (int i = 0; i < ktrtor::maxntt; ++i)
      ibuf[i] = ktrtor::tnx[i];
    tnx_vec.copyin(ibuf.data(), ktrtor::maxntt);
    for (int i = 0; i < ktrtor::maxntt; ++i)
      ibuf[i] = ktrtor::tny[i];
    tny_vec.copyin(ibuf.data(), ktrtor::maxntt);
    int count = ktrtor::maxtgrd * ktrtor::maxntt;
    ttx_vec.copyin(&ktrtor::ttx[0][0], count);
    tty_vec.copyin(&ktrtor::tty[0][0], count);
    count = ktrtor::maxtgrd2 * ktrtor::maxntt;
    tbf_vec.copyin(&ktrtor::tbf[0][0], count);
    tbx_vec.copyin(&ktrtor::tbx[0][0], count);
    tby_vec.copyin(&ktrtor::tby[0][0], count);
    tbxy_vec.copyin(&ktrtor::tbxy[0][0], count);

    ttorunit = torpot::ttorunit;

    // see also subroutine chkttor in etortor.f
    ibuf.resize(ntortor);
    for (int itortor = 0; itortor < ntortor; ++itortor) {
      int ib, ic, id;
      int i = tortor::itt[itortor * 3] - 1;
      int flag = tortor::itt[itortor * 3 + 2];
      if (flag == 1) {
        ib = bitor_::ibitor[5 * i + 1] - 1;
        ic = bitor_::ibitor[5 * i + 2] - 1;
        id = bitor_::ibitor[5 * i + 3] - 1;
      } else {
        ib = bitor_::ibitor[5 * i + 3] - 1;
        ic = bitor_::ibitor[5 * i + 2] - 1;
        id = bitor_::ibitor[5 * i + 1] - 1;
      }

      int ia = -1; // -1 is the default ia value for non-chiral center ic
      if (couple::n12[ic] == 4) {
        /*
         *  j    k
         *   \  /
         *    ic
         *   /  \
         * ib    id
         */
        int j = -1;
        int k;
        for (int i = 0; i < 4; ++i) {
          int m = couple::i12[ic][i] - 1;
          if (m != ib && m != id) {
            if (j == -1) // if j is uninitialized
              j = m;
            else // j has been uninitialized
              k = m;
          }
        }

        if (atoms::type[j] > atoms::type[k])
          ia = j;
        if (atoms::type[k] > atoms::type[j])
          ia = k;
        // else ia is uninitialized
        if (atomid::atomic[j] > atomid::atomic[k])
          ia = j;
        if (atomid::atomic[k] > atomid::atomic[j])
          ia = k;
        // else ia is unchanged
      }
      ibuf[itortor] = ia;
    }
    chkttor_ia_vec_.copyin(ibuf.data(), ntortor);
  }
}

extern void etortor_acc_impl_(int vers);
void etortor(int vers) { etortor_acc_impl_(vers); }
TINKER_NAMESPACE_END
