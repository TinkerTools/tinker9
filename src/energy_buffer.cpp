#include "energy_buffer.h"
#include "box.h"
#include "e_vdw.h"
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


void buffer_allocate(energy_buffer* pe, virial_buffer* pv)
{
   if (rc_flag & calc::analyz) {
      auto len = buffer_size();
      device_array::allocate(len, pe, pv);
      energy_buffers.push_back(*pe);
      virial_buffers.push_back(*pv);
   } else {
      *pe = esum_buf;
      *pv = vir_buf;
   }
}


void buffer_deallocate(energy_buffer e, virial_buffer v)
{
   if (rc_flag & calc::analyz)
      device_array::deallocate(e, v);
}


void buffer_allocate(count_buffer* pc, energy_buffer* pe, virial_buffer* pv)
{
   if (rc_flag & calc::analyz) {
      auto len = buffer_size();
      device_array::allocate(len, pc, pe, pv);
      count_buffers.push_back(*pc);
      energy_buffers.push_back(*pe);
      virial_buffers.push_back(*pv);
   } else {
      *pc = nullptr;
      *pe = esum_buf;
      *pv = vir_buf;
   }
}


void buffer_deallocate(count_buffer c, energy_buffer e, virial_buffer v)
{
   if (rc_flag & calc::analyz)
      device_array::deallocate(c, e, v);
}


int get_count(const count_buffer ne)
{
   int c = parallel::reduce_sum(ne, buffer_size(), WAIT_NEW_Q);
   return c;
}


real get_energy(const energy_buffer e)
{
   auto b = parallel::reduce_sum(e, buffer_size(), WAIT_NEW_Q);
   real real_out = energy_buffer_traits::cast(b);


   // vdw long-range correction
   // check != 0 for non-PBC
   if (e == ev && elrc_vol != 0) {
      real_out += elrc_vol / volbox();
   }


   return real_out;
}


void get_virial(real (&v1)[virial_buffer_traits::n], const virial_buffer v)
{
   virial_buffer_traits::type b[virial_buffer_traits::n];
   parallel::reduce_sum2(b, v, buffer_size(), WAIT_NEW_Q);
   for (size_t i = 0; i < virial_buffer_traits::n; ++i)
      v1[i] = virial_buffer_traits::cast(b[i]);
}


void get_virial(real (&v_out)[9], const virial_buffer v)
{
   real v1[virial_buffer_traits::n];
   get_virial(v1, v);
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
      real term = vlrc_vol / volbox();
      v_out[0] += term; // xx
      v_out[4] += term; // yy
      v_out[8] += term; // zz
   }
}


std::vector<count_buffer> count_buffers;
std::vector<energy_buffer> energy_buffers;
std::vector<virial_buffer> virial_buffers;
TINKER_NAMESPACE_END
