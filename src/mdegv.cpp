#include "mdegv.h"
#include "mdcalc.h"
#include "mdpq.h"


namespace tinker {
energy_prec esum, eksum, ekin[3][3];
energy_buffer eng_buf;
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
         darray::zero(PROCEED_NEW_Q, buffer_size(), u);
      }
   }

   if (vers & calc::energy) {
      esum = 0;
      size_t size = energy_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         energy_buffer u = energy_buffers[i];
         darray::zero(PROCEED_NEW_Q, buffer_size(), u);
      }
   }

   if (vers & calc::grad) {
      size_t size = x_grads.size();
      for (size_t i = 0; i < size; ++i) {
         darray::zero(PROCEED_NEW_Q, n, x_grads[i], y_grads[i], z_grads[i]);
      }
   }

   if (vers & calc::virial) {
      for (int i = 0; i < 9; ++i) {
         vir[i] = 0;
      }
      size_t size = virial_buffers.size();
      for (size_t i = 0; i < size; ++i) {
         virial_buffer u = virial_buffers[i];
         darray::zero(PROCEED_NEW_Q, buffer_size(), u);
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


void scale_gradient(double scale, grad_prec* g0x, grad_prec* g0y,
                    grad_prec* g0z)
{
   if (scale == 1)
      return;
   else if (scale == 0) {
      zero_gradient(PROCEED_NEW_Q, n, g0x, g0y, g0z);
   } else
      scale_gradient_acc(scale, g0x, g0y, g0z);
}


void sum_gradient(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
                  const grad_prec* g1x, const grad_prec* g1y,
                  const grad_prec* g1z)
{
   sum_gradient_acc(g0x, g0y, g0z, g1x, g1y, g1z);
}


void sum_gradient(double s, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
                  const grad_prec* g1x, const grad_prec* g1y,
                  const grad_prec* g1z)
{
   sum_gradient_acc(s, g0x, g0y, g0z, g1x, g1y, g1z);
}


//====================================================================//


void sum_energy(int vers)
{
   if (vers & calc::energy) {
      for (size_t i = 0; i < energy_buffers.size(); ++i) {
         energy_buffer u = energy_buffers[i];
         energy_prec e = energy_reduce(u);
         energy_prec* eptr = get_energy_reduce_dst(u);
         *eptr = e;
         if (eptr != &esum)
            esum += e;
      }
   }

   if (vers & calc::virial) {
      for (size_t i = 0; i < virial_buffers.size(); ++i) {
         virial_buffer u = virial_buffers[i];
         virial_prec v[9];
         virial_reduce(v, u);
         for (int iv = 0; iv < 9; ++iv)
            vir[iv] += v[iv];
      }
   }

   if (vers & calc::grad) {
      size_t ngrad = x_grads.size();
      for (size_t i = 1; i < ngrad; ++i) {
         sum_gradient(gx, gy, gz, x_grads[i], y_grads[i], z_grads[i]);
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
         darray::copyout(PROCEED_NEW_Q, n, hgx.data(), gx_src);
         darray::copyout(PROCEED_NEW_Q, n, hgy.data(), gy_src);
         darray::copyout(WAIT_NEW_Q, n, hgz.data(), gz_src);
         for (int i = 0; i < n; ++i) {
            grdx[i] = to_flt_host<double>(hgx[i]);
            grdy[i] = to_flt_host<double>(hgy[i]);
            grdz[i] = to_flt_host<double>(hgz[i]);
         }
#else
         if (sizeof(grad_prec) < sizeof(double)) {
            std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
            darray::copyout(PROCEED_NEW_Q, n, hgx.data(), gx_src);
            darray::copyout(PROCEED_NEW_Q, n, hgy.data(), gy_src);
            darray::copyout(WAIT_NEW_Q, n, hgz.data(), gz_src);
            for (int i = 0; i < n; ++i) {
               grdx[i] = hgx[i];
               grdy[i] = hgy[i];
               grdz[i] = hgz[i];
            }
         } else {
            darray::copyout(PROCEED_NEW_Q, n, grdx, (double*)gx_src);
            darray::copyout(PROCEED_NEW_Q, n, grdy, (double*)gy_src);
            darray::copyout(WAIT_NEW_Q, n, grdz, (double*)gz_src);
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


void egv_data(rc_op op)
{
   if (rc_flag & calc::energy) {
      if (op & rc_dealloc) {
         if (rc_flag & calc::analyz) {
         } else {
            darray::deallocate(eng_buf);
         }
      }


      if (op & rc_alloc) {
         count_buffers.clear();
         energy_buffers.clear();
         clear_energy_reduce_dst();
         if (rc_flag & calc::analyz) {
            eng_buf = nullptr;
         } else {
            auto sz = buffer_size();
            darray::allocate(sz, &eng_buf);
            energy_buffers.push_back(eng_buf);
            set_energy_reduce_dst(eng_buf, &esum);
         }
      }
   }


   if (rc_flag & calc::virial) {
      if (op & rc_dealloc) {
         if (rc_flag & calc::analyz) {
         } else {
            darray::deallocate(vir_buf);
         }
      }


      if (op & rc_alloc) {
         virial_buffers.clear();
         if (rc_flag & calc::analyz) {
            vir_buf = nullptr;
         } else {
            auto sz = buffer_size();
            darray::allocate(sz, &vir_buf);
            virial_buffers.push_back(vir_buf);
         }
      }
   }


   if (rc_flag & calc::grad) {
      if (op & rc_dealloc) {
         // Always deallocate gx, gy, gz.
         // Other gradients are deallocated elsewhere.
         darray::deallocate(gx, gy, gz);
      }


      if (op & rc_alloc) {
         // Always allocate gx, gy, gz as the first gradient.
         x_grads.clear();
         y_grads.clear();
         z_grads.clear();
         darray::allocate(n, &gx, &gy, &gz);
         x_grads.push_back(gx);
         y_grads.push_back(gy);
         z_grads.push_back(gz);
      }
   }
}
}
