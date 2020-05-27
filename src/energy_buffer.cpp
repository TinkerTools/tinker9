#include "energy_buffer.h"
#include "box.h"
#include "evdw.h"
#include "mathfunc.h"
#include "md.h"
#include <cassert>

namespace tinker {
size_t buffer_size()
{
   size_t nhint = n;
   size_t max_bytes = 2 * 1024 * 1024ull; // 2 MB for int
   if (nhint <= 16384)
      max_bytes /= 2; // 1 MB
   if (nhint <= 8192)
      max_bytes /= 2; // 512 KB
   if (nhint <= 4096)
      max_bytes /= 2; // 256 KB
   if (nhint <= 2048)
      max_bytes /= 2; // 128 KB
   if (nhint <= 1024)
      max_bytes /= 2; // 64 KB
   if (nhint <= 512)
      max_bytes /= 2; // 32 KB 8192 words

   assert(is_pow2(max_bytes) && "new_size must be power of 2");
   return max_bytes;
}


void buffer_allocate(int flag, energy_buffer* pe, virial_buffer* pv)
{
   if (flag & calc::analyz) {
      auto len = buffer_size();
      if (flag & calc::energy) {
         darray::allocate(len, pe);
      }
      if (flag & calc::virial) {
         darray::allocate(len, pv);
      }
   } else {
      if (flag & calc::energy) {
         *pe = eng_buf;
      }
      if (flag & calc::virial) {
         *pv = vir_buf;
      }
   }
}


void buffer_deallocate(int flag, energy_buffer e, virial_buffer v)
{
   if (flag & calc::analyz) {
      if (flag & calc::energy) {
         darray::deallocate(e);
      }
      if (flag & calc::virial) {
         darray::deallocate(v);
      }
   }
}


void buffer_allocate(int flag, count_buffer* pc)
{
   if (flag & calc::analyz) {
      auto len = buffer_size();
      darray::allocate(len, pc);
   } else {
      *pc = nullptr;
   }
}


void buffer_deallocate(int flag, count_buffer c)
{
   if (flag & calc::analyz) {
      darray::deallocate(c);
   }
}


void buffer_allocate(int flag, grad_prec** px, grad_prec** py, grad_prec** pz)
{
   if (flag & calc::grad) {
      if (flag & calc::analyz) {
         darray::allocate(n, px, py, pz);
      } else {
         *px = gx;
         *py = gy;
         *pz = gz;
      }
   }
}


void buffer_deallocate(int flag, grad_prec* gx, grad_prec* gy, grad_prec* gz)
{
   if (flag & calc::grad) {
      if (flag & calc::analyz) {
         darray::deallocate(gx, gy, gz);
      }
   }
}

int count_reduce(const count_buffer ne)
{
   int c = parallel::reduce_sum(ne, buffer_size(), WAIT_NEW_Q);
   return c;
}


energy_prec energy_reduce(const energy_buffer e)
{
   auto b = parallel::reduce_sum(e, buffer_size(), WAIT_NEW_Q);
   energy_prec real_out = to_flt_host<energy_prec>(b);
   return real_out;
}


void virial_reduce(virial_prec (&v1)[virial_buffer_traits::N],
                   const virial_buffer v)
{
   virial_buffer_traits::type b[virial_buffer_traits::N];
   parallel::reduce_sum2(b, v, buffer_size(), WAIT_NEW_Q);
   for (size_t i = 0; i < virial_buffer_traits::N; ++i)
      v1[i] = to_flt_host<virial_prec>(b[i]);
}


void virial_reduce(virial_prec (&v_out)[9], const virial_buffer v)
{
   virial_prec v1[virial_buffer_traits::N];
   virial_reduce(v1, v);
   // xx yx zx yy zy zz
   //  0  1  2  3  4  5
   v_out[0] = v1[0]; // xx
   v_out[1] = v1[1]; // xy
   v_out[2] = v1[2]; // xz
   v_out[3] = v1[1]; // yx
   v_out[4] = v1[3]; // yy
   v_out[5] = v1[4]; // yz
   v_out[6] = v1[2]; // zx
   v_out[7] = v1[4]; // zy
   v_out[8] = v1[5]; // zz
}
}
