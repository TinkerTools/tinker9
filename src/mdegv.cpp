#include "mdegv.h"
#include "host_zero.h"
#include "mdcalc.h"
#include "mdpq.h"
#include "mod.energi.h"


namespace tinker {
energy_prec eksum, ekin[3][3];


//====================================================================//


void zero_egv(int vers)
{
   if (vers & calc::analyz) {
      darray::zero(PROCEED_NEW_Q, buffer_size(), TINKER_COUNT_BUFFERS);
   }

   if (vers & calc::energy) {
      host_zero(TINKER_ENERGY_VARIABLES);
      darray::zero(PROCEED_NEW_Q, buffer_size(), TINKER_ENERGY_BUFFERS);
   }

   if (vers & calc::virial) {
      host_zero(TINKER_VIRIAL_TENSORS);
      darray::zero(PROCEED_NEW_Q, buffer_size(), TINKER_VIRIAL_BUFFERS);
   }

   if (vers & calc::grad) {
      darray::zero(PROCEED_NEW_Q, n, TINKER_GRADIENTS);
   }
}


void zero_egv()
{
   zero_egv(rc_flag);
}


//====================================================================//


void scale_gradient(double scale, grad_prec* g0x, grad_prec* g0y,
                    grad_prec* g0z)
{
   if (scale == 1)
      return;
   else if (scale == 0) {
      darray::zero(PROCEED_NEW_Q, n, g0x, g0y, g0z);
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
         host_zero(TINKER_COUNT_BUFFERS);
         host_zero(TINKER_ENERGY_BUFFERS);
         if (rc_flag & calc::analyz) {
            eng_buf = nullptr;
         } else {
            auto sz = buffer_size();
            darray::allocate(sz, &eng_buf);
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
         host_zero(TINKER_VIRIAL_BUFFERS);
         if (rc_flag & calc::analyz) {
            vir_buf = nullptr;
         } else {
            auto sz = buffer_size();
            darray::allocate(sz, &vir_buf);
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
         host_zero(TINKER_GRADIENTS);
         darray::allocate(n, &gx, &gy, &gz);
      }
   }
}
}
