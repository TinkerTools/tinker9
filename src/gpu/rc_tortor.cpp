#include "gpu/decl_mdstate.h"
#include "gpu/e_tortor.h"
#include "gpu/rc.h"
#include "util_potent.h"

TINKER_NAMESPACE_BEGIN
// module bitor
int nbitor;
int (*ibitor)[5];

// module tortor
int ntortor;
int (*itt)[3];

// module ktrtor
int *tnx, *tny;
real (*ttx)[ktrtor::maxtgrd];
real (*tty)[ktrtor::maxtgrd];
real (*tbf)[ktrtor::maxtgrd2];
real (*tbx)[ktrtor::maxtgrd2];
real (*tby)[ktrtor::maxtgrd2];
real (*tbxy)[ktrtor::maxtgrd2];

real ttorunit;

int* chkttor_ia_;

real* ett;
real* vir_ett;

void etortor_data(rc_t rc) {
  if (!use_potent(tortor_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(ibitor));

    check_cudart(cudaFree(itt));

    check_cudart(cudaFree(tnx));
    check_cudart(cudaFree(tny));
    check_cudart(cudaFree(ttx));
    check_cudart(cudaFree(tty));
    check_cudart(cudaFree(tbf));
    check_cudart(cudaFree(tbx));
    check_cudart(cudaFree(tby));
    check_cudart(cudaFree(tbxy));

    check_cudart(cudaFree(chkttor_ia_));

    free_ev(ett, vir_ett);
  }

  if (rc & rc_alloc) {
    nbitor = bitor_::nbitor;
    check_cudart(cudaMalloc(&ibitor, sizeof(int) * 5 * nbitor));

    check_cudart(cudaMalloc(&itt, sizeof(int) * 3 * nbitor));

    const size_t rs = sizeof(real);
    check_cudart(cudaMalloc(&tnx, sizeof(int) * ktrtor::maxntt));
    check_cudart(cudaMalloc(&tny, sizeof(int) * ktrtor::maxntt));
    int count = ktrtor::maxtgrd * ktrtor::maxntt;
    check_cudart(cudaMalloc(&ttx, rs * count));
    check_cudart(cudaMalloc(&tty, rs * count));
    count = ktrtor::maxtgrd2 * ktrtor::maxntt;
    check_cudart(cudaMalloc(&tbf, rs * count));
    check_cudart(cudaMalloc(&tbx, rs * count));
    check_cudart(cudaMalloc(&tby, rs * count));
    check_cudart(cudaMalloc(&tbxy, rs * count));

    ntortor = count_bonded_term(tortor_term);
    check_cudart(cudaMalloc(&chkttor_ia_, sizeof(int) * ntortor));

    alloc_ev(&ett, &vir_ett);
  }

  if (rc & rc_copyin) {
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
