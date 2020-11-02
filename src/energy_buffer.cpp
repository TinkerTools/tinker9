#include "tool/energy_buffer.h"
#include "box.h"
#include "evdw.h"
#include "mathfunc.h"
#include "md.h"
#include <cassert>

namespace tinker {
size_t buffer_size()
{
   size_t bsize = nelem_buffer;
   assert(is_pow2(bsize) && "buffer size must be power of 2.");
   return bsize;
}


void buffer_allocate(int flag, energy_buffer* pe, virial_buffer* pv,
                     grad_prec** px, grad_prec** py, grad_prec** pz)
{
   assert(flag & calc::analyz);
   size_t bsize = buffer_size();
   if (flag & calc::energy)
      darray::allocate(bsize, pe);
   if (flag & calc::virial)
      darray::allocate(bsize, pv);
   if (flag & calc::grad)
      darray::allocate(n, px, py, pz);
}


void buffer_deallocate(int flag, energy_buffer e, virial_buffer v,
                       grad_prec* gx, grad_prec* gy, grad_prec* gz)
{
   assert(flag & calc::analyz);
   if (flag & calc::energy)
      darray::deallocate(e);
   if (flag & calc::virial)
      darray::deallocate(v);
   if (flag & calc::grad)
      darray::deallocate(gx, gy, gz);
}


void buffer_allocate(int flag, count_buffer* pc)
{
   assert(flag & calc::analyz);
   size_t bsize = buffer_size();
   darray::allocate(bsize, pc);
}


void buffer_deallocate(int flag, count_buffer c)
{
   assert(flag & calc::analyz);
   darray::deallocate(c);
}
int count_reduce(const count_buffer ne)
{
   int c = parallel::reduce_sum(ne, buffer_size(), asyncq);
   return c;
}


energy_prec energy_reduce(const energy_buffer e)
{
   auto b = parallel::reduce_sum(e, buffer_size(), asyncq);
   energy_prec real_out = to_flt_host<energy_prec>(b);
   return real_out;
}


void virial_reduce(virial_prec (&v1)[virial_buffer_traits::N],
                   const virial_buffer v)
{
   virial_buffer_traits::type b[virial_buffer_traits::N];
   parallel::reduce_sum2(b, v, buffer_size(), asyncq);
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
