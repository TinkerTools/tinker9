#include "couple.h"
#include "md.h"
#include "tool/darray.h"
#include <tinker/detail/couple.hh>
#include <tinker/detail/sizes.hh>
#include <vector>


namespace tinker {
static_assert(couple_maxn12 >= sizes::maxval, "");


int (*couple_i12)[couple_maxn12];
int* couple_n12;


void couple_data(rc_op op)
{
   if (op & rc_dealloc) {
      darray::deallocate(couple_i12, couple_n12);
   }


   if (op & rc_alloc) {
      darray::allocate(n, &couple_i12, &couple_n12);
   }


   if (op & rc_init) {
      std::vector<int> ibuf;
      ibuf.resize(couple_maxn12 * n);
      for (int i = 0; i < n; ++i) {
         int nn = couple::n12[i];
         int base = i * couple_maxn12;
         for (int j = 0; j < nn; ++j) {
            int k = couple::i12[i][j];
            ibuf[base + j] = k - 1;
         }
         for (int j = nn; j < couple_maxn12; ++j) {
            ibuf[base + j] = -1;
         }
      }
      darray::copyin(g::q0, n, couple_i12, ibuf.data());
      darray::copyin(g::q0, n, couple_n12, couple::n12);
      wait_for(g::q0);
   }
}
}
