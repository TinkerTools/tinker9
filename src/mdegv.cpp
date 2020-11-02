#include "mdegv.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdpq.h"
#include "tool/device_zero.h"
#include "tool/host_zero.h"


namespace tinker {


//====================================================================//


void zero_egv(int vers)
{
   size_t bsize = buffer_size();
   if (vers & calc::energy) {
      host_zero(esum, energy_valence, energy_vdw, energy_elec);
      zero3_async(bsize, eng_buf, eng_buf_vdw, eng_buf_elec);
   }

   if (vers & calc::virial) {
      host_zero(vir, virial_valence, virial_vdw, virial_elec);
      zero3_async(bsize, vir_buf, vir_buf_vdw, vir_buf_elec);
   }

   if (vers & calc::grad) {
      zero9_async(n, gx, gy, gz, gx_vdw, gy_vdw, gz_vdw, gx_elec, gy_elec,
                  gz_elec);
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
      darray::zero(async_queue, n, g0x, g0y, g0z);
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
                   const grad_prec* gz_src, int queue)
{
   if (vers & calc::grad) {
      if (grdx && grdy && grdz) {
#if TINKER_DETERMINISTIC_FORCE
         std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
         darray::copyout(queue, n, hgx.data(), gx_src);
         darray::copyout(queue, n, hgy.data(), gy_src);
         darray::copyout(queue, n, hgz.data(), gz_src);
         wait_for(queue);
         for (int i = 0; i < n; ++i) {
            grdx[i] = to_flt_host<double>(hgx[i]);
            grdy[i] = to_flt_host<double>(hgy[i]);
            grdz[i] = to_flt_host<double>(hgz[i]);
         }
#else
         if (sizeof(grad_prec) < sizeof(double)) {
            std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
            darray::copyout(queue, n, hgx.data(), gx_src);
            darray::copyout(queue, n, hgy.data(), gy_src);
            darray::copyout(queue, n, hgz.data(), gz_src);
            wait_for(queue);
            for (int i = 0; i < n; ++i) {
               grdx[i] = hgx[i];
               grdy[i] = hgy[i];
               grdz[i] = hgz[i];
            }
         } else {
            darray::copyout(queue, n, grdx, (double*)gx_src);
            darray::copyout(queue, n, grdy, (double*)gy_src);
            darray::copyout(queue, n, grdz, (double*)gz_src);
            wait_for(queue);
         }
#endif
      }
   }
}


void copy_gradient(int vers, double* grdx, double* grdy, double* grdz,
                   const grad_prec* gx_src, const grad_prec* gy_src,
                   const grad_prec* gz_src)
{
   copy_gradient(vers, grdx, grdy, grdz, gx_src, gy_src, gz_src, async_queue);
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
   bool rc_a = rc_flag & calc::analyz;


   if (rc_flag & calc::energy) {
      if (op & rc_dealloc) {
         if (!rc_a) {
            darray::deallocate(eng_buf);
            if (use_energi_vdw())
               darray::deallocate(eng_buf_vdw);
            if (use_energi_elec())
               darray::deallocate(eng_buf_elec);
         }
      }


      if (op & rc_alloc) {
         host_zero(eng_buf, eng_buf_vdw, eng_buf_elec);
         if (!rc_a) {
            auto sz = buffer_size();
            darray::allocate(sz, &eng_buf);
            if (use_energi_vdw())
               darray::allocate(sz, &eng_buf_vdw);
            if (use_energi_elec())
               darray::allocate(sz, &eng_buf_elec);
         }
      }
   }


   if (rc_flag & calc::virial) {
      if (op & rc_dealloc) {
         if (!rc_a) {
            darray::deallocate(vir_buf);
            if (use_energi_vdw())
               darray::deallocate(vir_buf_vdw);
            if (use_energi_elec())
               darray::deallocate(vir_buf_elec);
         }
      }


      if (op & rc_alloc) {
         host_zero(vir_buf, vir_buf_vdw, vir_buf_elec);
         if (!rc_a) {
            auto sz = buffer_size();
            darray::allocate(sz, &vir_buf);
            if (use_energi_vdw())
               darray::allocate(sz, &vir_buf_vdw);
            if (use_energi_elec())
               darray::allocate(sz, &vir_buf_elec);
         }
      }
   }


   if (rc_flag & calc::grad) {
      if (op & rc_dealloc) {
         darray::deallocate(gx, gy, gz);
         if (use_energi_vdw())
            darray::deallocate(gx_vdw, gy_vdw, gz_vdw);
         if (use_energi_elec())
            darray::deallocate(gx_elec, gy_elec, gz_elec);
      }


      if (op & rc_alloc) {
         host_zero(gx, gy, gz, gx_vdw, gy_vdw, gz_vdw, gx_elec, gy_elec,
                   gz_elec);
         darray::allocate(n, &gx, &gy, &gz);
         if (use_energi_vdw())
            darray::allocate(n, &gx_vdw, &gy_vdw, &gz_vdw);
         if (use_energi_elec())
            darray::allocate(n, &gx_elec, &gy_elec, &gz_elec);
      }
   }
}
}
