#include "e_tortor.h"
#include "array.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/atomid.hh>
#include <ext/tinker/detail/atoms.hh>
#include <ext/tinker/detail/bitor.hh>
#include <ext/tinker/detail/couple.hh>
#include <ext/tinker/detail/torpot.hh>
#include <ext/tinker/detail/tortor.hh>

TINKER_NAMESPACE_BEGIN
void etortor_data(rc_op op) {
  if (!use_potent(tortor_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(ibitor);

    dealloc_bytes(itt);

    dealloc_bytes(tnx);
    dealloc_bytes(tny);
    dealloc_bytes(ttx);
    dealloc_bytes(tty);
    dealloc_bytes(tbf);
    dealloc_bytes(tbx);
    dealloc_bytes(tby);
    dealloc_bytes(tbxy);

    dealloc_bytes(chkttor_ia_);

    dealloc_ev(ett, vir_ett);
  }

  if (op & rc_alloc) {
    nbitor = bitor_::nbitor;
    alloc_bytes(&ibitor, sizeof(int) * 5 * nbitor);

    alloc_bytes(&itt, sizeof(int) * 3 * nbitor);

    const size_t rs = sizeof(real);
    alloc_bytes(&tnx, sizeof(int) * ktrtor::maxntt);
    alloc_bytes(&tny, sizeof(int) * ktrtor::maxntt);
    int count = ktrtor::maxtgrd * ktrtor::maxntt;
    alloc_bytes(&ttx, rs * count);
    alloc_bytes(&tty, rs * count);
    count = ktrtor::maxtgrd2 * ktrtor::maxntt;
    alloc_bytes(&tbf, rs * count);
    alloc_bytes(&tbx, rs * count);
    alloc_bytes(&tby, rs * count);
    alloc_bytes(&tbxy, rs * count);

    ntortor = count_bonded_term(tortor_term);
    alloc_bytes(&chkttor_ia_, sizeof(int) * ntortor);

    alloc_ev(&ett, &vir_ett);
  }

  if (op & rc_init) {
    std::vector<int> ibuf;

    ibuf.resize(5 * nbitor);
    for (int i = 0; i < 5 * nbitor; ++i)
      ibuf[i] = bitor_::ibitor[i] - 1;
    copyin_array(&ibitor[0][0], ibuf.data(), 5 * nbitor);

    ibuf.resize(3 * nbitor);
    for (int i = 0; i < 3 * nbitor; ++i)
      ibuf[i] = tortor::itt[i] - 1;
    copyin_array(&itt[0][0], ibuf.data(), 3 * nbitor);

    ibuf.resize(ktrtor::maxntt);
    for (int i = 0; i < ktrtor::maxntt; ++i)
      ibuf[i] = ktrtor::tnx[i];
    copyin_array(tnx, ibuf.data(), ktrtor::maxntt);
    for (int i = 0; i < ktrtor::maxntt; ++i)
      ibuf[i] = ktrtor::tny[i];
    copyin_array(tny, ibuf.data(), ktrtor::maxntt);
    int count = ktrtor::maxtgrd * ktrtor::maxntt;
    copyin_array(&ttx[0][0], &ktrtor::ttx[0][0], count);
    copyin_array(&tty[0][0], &ktrtor::tty[0][0], count);
    count = ktrtor::maxtgrd2 * ktrtor::maxntt;
    copyin_array(&tbf[0][0], &ktrtor::tbf[0][0], count);
    copyin_array(&tbx[0][0], &ktrtor::tbx[0][0], count);
    copyin_array(&tby[0][0], &ktrtor::tby[0][0], count);
    copyin_array(&tbxy[0][0], &ktrtor::tbxy[0][0], count);

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
    copyin_array(chkttor_ia_, ibuf.data(), ntortor);
  }
}

extern void etortor_acc_impl_(int vers);
void etortor(int vers) { etortor_acc_impl_(vers); }
TINKER_NAMESPACE_END
