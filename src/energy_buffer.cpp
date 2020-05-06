#include "energy_buffer.h"
#include "box.h"
#include "evdw.h"
#include "mathfunc.h"
#include "md.h"
#include <cassert>

TINKER_NAMESPACE_BEGIN
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


void buffer_allocate(int flag, energy_buffer* pe, grad_prec** px,
                     grad_prec** py, grad_prec** pz, virial_buffer* pv,
                     energy_prec* eptr)
{
   if (flag & calc::analyz) {
      auto len = buffer_size();
      if (flag & calc::energy) {
         darray::allocate(len, pe);
         energy_buffers.push_back(*pe);
         set_energy_reduce_dst(*pe, eptr);
      }
      if (flag & calc::grad) {
         darray::allocate(n, px, py, pz);
         x_grads.push_back(*px);
         y_grads.push_back(*py);
         z_grads.push_back(*pz);
      }
      if (flag & calc::virial) {
         darray::allocate(len, pv);
         virial_buffers.push_back(*pv);
      }
   } else {
      if (flag & calc::energy) {
         *pe = eng_buf;
      }
      if (flag & calc::grad) {
         *px = gx;
         *py = gy;
         *pz = gz;
      }
      if (flag & calc::virial) {
         *pv = vir_buf;
      }
   }
}


void buffer_deallocate(int flag, energy_buffer e, grad_prec* gx, grad_prec* gy,
                       grad_prec* gz, virial_buffer v)
{
   if (flag & calc::analyz) {
      if (flag & calc::energy) {
         darray::deallocate(e);
      }
      if (flag & calc::grad) {
         darray::deallocate(gx, gy, gz);
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
      count_buffers.push_back(*pc);
   } else {
      *pc = nullptr;
   }
}


void buffer_deallocate(int flag, count_buffer c)
{
   if (flag & calc::analyz)
      darray::deallocate(c);
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


   // vdw long-range correction
   // check != 0 for non-PBC
   if (e == ev && elrc_vol != 0) {
      real_out += elrc_vol / volbox();
   }


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


   // vdw long-range correction
   // check != 0 for non-PBC
   if (v == vir_ev && vlrc_vol != 0) {
      virial_prec term = vlrc_vol / volbox();
      v_out[0] += term; // xx
      v_out[4] += term; // yy
      v_out[8] += term; // zz
   }
}


namespace {
std::map<energy_buffer, energy_prec*> edst;
}


void set_energy_reduce_dst(energy_buffer b, energy_prec* d)
{
   edst[b] = d;
}


energy_prec* get_energy_reduce_dst(energy_buffer b)
{
   return edst.at(b);
}


void clear_energy_reduce_dst()
{
   edst.clear();
}


std::vector<count_buffer> count_buffers;
std::vector<energy_buffer> energy_buffers;
std::vector<virial_buffer> virial_buffers;
std::vector<grad_prec*> x_grads, y_grads, z_grads;
TINKER_NAMESPACE_END
