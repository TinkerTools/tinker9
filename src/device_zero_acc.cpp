#include "tool/device_zero_acc.h"


namespace tinker {
template <class T>
void zero3_async_acc(int nelem, T* a1, T* a2, T* a3)
{
   #pragma acc parallel loop async deviceptr(a1,a2,a3)
   for (int i = 0; i < nelem; ++i) {
      a1[i] = 0, a2[i] = 0, a3[i] = 0;
   }
}
template void zero3_async_acc(int, fixed*, fixed*, fixed*);
template void zero3_async_acc(int, real*, real*, real*);


template <class T, int N>
void zero3_async_acc(int nelem, T (*b1)[N], T (*b2)[N], T (*b3)[N])
{
   T* a1 = (T*)b1;
   T* a2 = (T*)b2;
   T* a3 = (T*)b3;
   #pragma acc parallel loop async deviceptr(a1,a2,a3)
   for (int i = 0; i < N * nelem; ++i) {
      a1[i] = 0, a2[i] = 0, a3[i] = 0;
   }
}
template void zero3_async_acc(int, fixed (*)[8], fixed (*)[8], fixed (*)[8]);
template void zero3_async_acc(int, real (*)[8], real (*)[8], real (*)[8]);


template <class T>
void zero9_async_acc(int nelem, T* a1, T* a2, T* a3, T* a4, T* a5, T* a6, T* a7,
                     T* a8, T* a9)
{
   #pragma acc parallel loop async deviceptr(a1,a2,a3,a4,a5,a6,a7,a8,a9)
   for (int i = 0; i < nelem; ++i) {
      a1[i] = 0, a2[i] = 0, a3[i] = 0, a4[i] = 0, a5[i] = 0;
      a6[i] = 0, a7[i] = 0, a8[i] = 0, a9[i] = 0;
   }
}
template void zero9_async_acc(int, fixed*, fixed*, fixed*, fixed*, fixed*,
                              fixed*, fixed*, fixed*, fixed*);
template void zero9_async_acc(int, real*, real*, real*, real*, real*, real*,
                              real*, real*, real*);
}
