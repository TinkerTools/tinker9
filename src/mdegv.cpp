#include "ff/energy.h"
#include "math/zero.h"

namespace tinker {
void zero_egv(int vers)
{
   size_t bsize = bufferSize();
   if (vers & calc::energy) {
      zeroOnHost(esum, energy_valence, energy_vdw, energy_elec);
      zeroOnDevice3Async(bsize, eng_buf, eng_buf_vdw, eng_buf_elec);
   }

   if (vers & calc::virial) {
      zeroOnHost(vir, virial_valence, virial_vdw, virial_elec);
      zeroOnDevice3Async(bsize, vir_buf, vir_buf_vdw, vir_buf_elec);
   }

   if (vers & calc::grad) {
      zeroOnDevice9Async(n, gx, gy, gz, gx_vdw, gy_vdw, gz_vdw, gx_elec, gy_elec, gz_elec);
   }
}

void zero_egv()
{
   zero_egv(rc_flag);
}

//====================================================================//

void scale_gradient(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z)
{
   if (scale == 1)
      return;
   else if (scale == 0) {
      darray::zero(g::q0, n, g0x, g0y, g0z);
   } else
      scale_gradient_acc(scale, g0x, g0y, g0z);
}

void sum_gradient(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z, const grad_prec* g1x,
   const grad_prec* g1y, const grad_prec* g1z)
{
   sum_gradient_acc(g0x, g0y, g0z, g1x, g1y, g1z);
}

void sum_gradient(double s, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z, const grad_prec* g1x,
   const grad_prec* g1y, const grad_prec* g1z)
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

void copy_gradient(int vers, double* grdx, double* grdy, double* grdz, const grad_prec* gx_src,
   const grad_prec* gy_src, const grad_prec* gz_src, int queue)
{
   if (vers & calc::grad) {
      if (grdx && grdy && grdz) {
#if TINKER_DETERMINISTIC_FORCE
         std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
         darray::copyout(queue, n, hgx.data(), gx_src);
         darray::copyout(queue, n, hgy.data(), gy_src);
         darray::copyout(queue, n, hgz.data(), gz_src);
         waitFor(queue);
         for (int i = 0; i < n; ++i) {
            grdx[i] = toFloat<double>(hgx[i]);
            grdy[i] = toFloat<double>(hgy[i]);
            grdz[i] = toFloat<double>(hgz[i]);
         }
#else
         if (sizeof(grad_prec) < sizeof(double)) {
            std::vector<grad_prec> hgx(n), hgy(n), hgz(n);
            darray::copyout(queue, n, hgx.data(), gx_src);
            darray::copyout(queue, n, hgy.data(), gy_src);
            darray::copyout(queue, n, hgz.data(), gz_src);
            waitFor(queue);
            for (int i = 0; i < n; ++i) {
               grdx[i] = hgx[i];
               grdy[i] = hgy[i];
               grdz[i] = hgz[i];
            }
         } else {
            darray::copyout(queue, n, grdx, (double*)gx_src);
            darray::copyout(queue, n, grdy, (double*)gy_src);
            darray::copyout(queue, n, grdz, (double*)gz_src);
            waitFor(queue);
         }
#endif
      }
   }
}

void copy_gradient(int vers, double* grdx, double* grdy, double* grdz, const grad_prec* gx_src,
   const grad_prec* gy_src, const grad_prec* gz_src)
{
   copy_gradient(vers, grdx, grdy, grdz, gx_src, gy_src, gz_src, g::q0);
}

void copy_gradient(int vers, double* grdx, double* grdy, double* grdz)
{
   copy_gradient(vers, grdx, grdy, grdz, gx, gy, gz);
}
}
