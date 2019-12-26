#include "energy_buffer.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
static bool use_ev_()
{
   return rc_flag & (calc::analyz | calc::energy | calc::virial);
}

static void grad_data_(rc_op op)
{
   if (!(rc_flag & calc::grad))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(gx, gy, gz);
   }

   if (op & rc_alloc) {
      device_array::allocate(n, &gx, &gy, &gz);
   }

   // we can never assume whether or not deriv::desum was allocated, because it
   // was allocated inside subroutine gradient(...), which would be skipped in
   // subroutine mdinit() if a dyn file existed to restart a simulation.

   // if (op & rc_init) {
   // copy in deriv::sum to gx, gy, and gz
   // }
}

void zero_egv()
{
   zero_egv(rc_flag & calc::vmask);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
static void ev_data_(rc_op op)
{
   if (!use_ev_())
      return;

   if (op & rc_dealloc) {
      if (!(rc_flag & calc::analyz))
         device_array::deallocate(esum_buf, vir_buf);
   }

   if (op & rc_alloc) {
      count_buffers.clear();
      energy_buffers.clear();
      virial_buffers.clear();
      if (rc_flag & calc::analyz) {
         esum_buf = nullptr;
         vir_buf = nullptr;
      } else {
         device_array::allocate(buffer_size(), &esum_buf, &vir_buf);
         energy_buffers.push_back(esum_buf);
         virial_buffers.push_back(vir_buf);
      }
   }
}

void egv_data(rc_op op)
{
   rc_man ev42_{ev_data_, op};
   rc_man grad42_{grad_data_, op};
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void zero_egv(int vers)
{
   if (vers & calc::analyz) {
      size_t size = count_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         count_buffer u = count_buffers[i];
         device_array::zero(buffer_size(), u);
      }
   }

   if (vers & calc::energy) {
      size_t size = energy_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         energy_buffer u = energy_buffers[i];
         device_array::zero(buffer_size(), u);
      }
   }

   if (vers & calc::virial) {
      size_t size = virial_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         virial_buffer u = virial_buffers[i];
         device_array::zero(buffer_size(), u);
      }
   }

   if (vers & calc::grad) {
      zero_gradient(n, gx, gy, gz);
   }
}

void sum_energies(int vers)
{
   if (vers & calc::energy) {
      if (rc_flag & calc::analyz) {
         esum = 0;
         for (size_t i = 0; i < energy_buffers.size(); ++i) {
            energy_buffer u = energy_buffers[i];
            real e = get_energy(u);
            esum += e;
         }
      } else {
         esum = get_energy(esum_buf);
      }
   }

   if (vers & calc::virial) {
      if (rc_flag & calc::analyz) {
         for (int iv = 0; iv < 9; ++iv)
            vir[iv] = 0;
         for (size_t i = 0; i < virial_buffers.size(); ++i) {
            virial_buffer u = virial_buffers[i];
            real v[9];
            get_virial(v, u);
            for (int iv = 0; iv < 9; ++iv)
               vir[iv] += v[iv];
         }
      } else {
         get_virial(vir, vir_buf);
      }
   }
}
TINKER_NAMESPACE_END
