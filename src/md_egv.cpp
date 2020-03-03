#include "md_egv.h"
#include "md_calc.h"
#include "md_pq.h"


TINKER_NAMESPACE_BEGIN
energy_prec esum, eksum, ekin[3][3];
energy_buffer esum_buf;
grad_prec *gx, *gy, *gz;
virial_prec vir[9];
virial_buffer vir_buf;


//====================================================================//


void zero_egv(int vers)
{
   if (vers & calc::analyz) {
      size_t size = count_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         count_buffer u = count_buffers[i];
         device_array::zero(PROCEED_NEW_Q, buffer_size(), u);
      }
   }

   if (vers & calc::energy) {
      size_t size = energy_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         energy_buffer u = energy_buffers[i];
         device_array::zero(PROCEED_NEW_Q, buffer_size(), u);
      }
   }

   if (vers & calc::virial) {
      size_t size = virial_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         virial_buffer u = virial_buffers[i];
         device_array::zero(PROCEED_NEW_Q, buffer_size(), u);
      }
   }

   if (vers & calc::grad) {
      zero_gradient(PROCEED_NEW_Q, n, gx, gy, gz);
   }
}


void zero_egv()
{
   zero_egv(rc_flag);
}


//====================================================================//


void zero_gradient(DMFlag flag, size_t nelem, real* gx, real* gy, real* gz)
{
   zero_gradient_acc(flag, nelem, gx, gy, gz);
}


void zero_gradient(DMFlag flag, size_t nelem, fixed* gx, fixed* gy, fixed* gz)
{
   zero_gradient_acc(flag, nelem, gx, gy, gz);
}


//====================================================================//


void sum_energy(int vers)
{
   if (vers & calc::energy) {
      if (rc_flag & calc::analyz) {
         esum = 0;
         for (size_t i = 0; i < energy_buffers.size(); ++i) {
            energy_buffer u = energy_buffers[i];
            energy_prec e = get_energy(u);
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
            virial_prec v[9];
            get_virial(v, u);
            for (int iv = 0; iv < 9; ++iv)
               vir[iv] += v[iv];
         }
      } else {
         get_virial(vir, vir_buf);
      }
   }
}


//====================================================================//


void copy_energy(int vers, energy_prec* eng)
{
   if (eng && vers & calc::energy && eng != &esum) {
      eng[0] = esum;
   }
}


void copy_gradient(int vers, double* grdx, double* grdy, double* grdz,
                   const grad_prec* gx_src, const grad_prec* gy_src,
                   const grad_prec* gz_src)
{
   if (vers & calc::grad) {
      if (grdx && grdy && grdz) {
#if TINKER_DETERMINISTIC_FORCE
         std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
         device_array::copyout(PROCEED_NEW_Q, n, hgx.data(), gx_src);
         device_array::copyout(PROCEED_NEW_Q, n, hgy.data(), gy_src);
         device_array::copyout(WAIT_NEW_Q, n, hgz.data(), gz_src);
         for (int i = 0; i < n; ++i) {
            grdx[i] = to_flt_host<double>(hgx[i]);
            grdy[i] = to_flt_host<double>(hgy[i]);
            grdz[i] = to_flt_host<double>(hgz[i]);
         }
#else
         if (sizeof(grad_prec) < sizeof(double)) {
            std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
            device_array::copyout(PROCEED_NEW_Q, n, hgx.data(), gx_src);
            device_array::copyout(PROCEED_NEW_Q, n, hgy.data(), gy_src);
            device_array::copyout(WAIT_NEW_Q, n, hgz.data(), gz_src);
            for (int i = 0; i < n; ++i) {
               grdx[i] = hgx[i];
               grdy[i] = hgy[i];
               grdz[i] = hgz[i];
            }
         } else {
            device_array::copyout(PROCEED_NEW_Q, n, grdx, (double*)gx_src);
            device_array::copyout(PROCEED_NEW_Q, n, grdy, (double*)gy_src);
            device_array::copyout(WAIT_NEW_Q, n, grdz, (double*)gz_src);
         }
#endif
      }
   }
}


void copy_gradient(int vers, double* grdx, double* grdy, double* grdz)
{
   copy_gradient(vers, grdx, grdy, grdz, gx, gy, gz);
}


void copy_virial(int vers, virial_prec* virial)
{
   if (virial && vers & calc::virial && virial != &vir[0]) {
      for (int i = 0; i < 9; ++i)
         virial[i] = vir[i];
   }
}


//====================================================================//


namespace {
bool use_ev()
{
   return rc_flag & (calc::analyz | calc::energy | calc::virial);
}


void grad_data(rc_op op)
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


void ev_data(rc_op op)
{
   if (!use_ev())
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
}


void egv_data(rc_op op)
{
   rc_man ev42_{ev_data, op};
   rc_man grad42_{grad_data, op};
}
TINKER_NAMESPACE_END
