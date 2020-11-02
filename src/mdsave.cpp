#include "mdsave.h"
#include "box.h"
#include "epolar.h"
#include "execq.h"
#include "potent.h"
#include "tinker_rt.h"
#include <condition_variable>
#include <future>
#include <mutex>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/deriv.hh>
#include <tinker/detail/moldyn.hh>
#include <tinker/detail/output.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/units.hh>


namespace tinker {
namespace {
std::mutex mtx_dup, mtx_write;
std::condition_variable cv_dup, cv_write;
bool idle_dup, idle_write;
std::future<void> fut_dup_then_write;


bool mdsave_use_uind()
{
   return output::uindsave && use_potent(polar_term);
}


ExecQ dup_stream;
real (*dup_buf_uind)[3];
energy_prec dup_buf_esum;
Box dup_buf_box;
pos_prec *dup_buf_x, *dup_buf_y, *dup_buf_z;
vel_prec *dup_buf_vx, *dup_buf_vy, *dup_buf_vz;
grad_prec *dup_buf_gx, *dup_buf_gy, *dup_buf_gz;


void mdsave_dup_then_write(int istep, time_prec dt)
{

   // duplicate

   dup_buf_esum = esum;
   get_default_box(dup_buf_box);
   darray::copy(asyncq, n, dup_buf_x, xpos);
   darray::copy(asyncq, n, dup_buf_y, ypos);
   darray::copy(asyncq, n, dup_buf_z, zpos);

   darray::copy(asyncq, n, dup_buf_vx, vx);
   darray::copy(asyncq, n, dup_buf_vy, vy);
   darray::copy(asyncq, n, dup_buf_vz, vz);

   darray::copy(asyncq, n, dup_buf_gx, gx);
   darray::copy(asyncq, n, dup_buf_gy, gy);
   darray::copy(asyncq, n, dup_buf_gz, gz);

   if (mdsave_use_uind()) {
      darray::copy(asyncq, 3 * n, &dup_buf_uind[0][0], &uind[0][0]);
   }

   // Record mdsave_begin_event when nonblk is available.
   // Stream (0) will wait until mdsave_begin_event is recorded.
   dup_stream.begin_copyout();

   mtx_dup.lock();
   idle_dup = true;
   cv_dup.notify_all();
   mtx_dup.unlock();

   // get gpu buffer and write to external files

   energy_prec epot = dup_buf_esum;
   set_tinker_box_module(dup_buf_box);
   if (sizeof(pos_prec) == sizeof(double)) {
      darray::copyout(syncq, n, atoms::x, dup_buf_x);
      darray::copyout(syncq, n, atoms::y, dup_buf_y);
      darray::copyout(syncq, n, atoms::z, dup_buf_z);
      wait_for(syncq);
   } else {
      std::vector<pos_prec> arrx(n), arry(n), arrz(n);
      darray::copyout(syncq, n, arrx.data(), dup_buf_x);
      darray::copyout(syncq, n, arry.data(), dup_buf_y);
      darray::copyout(syncq, n, arrz.data(), dup_buf_z);
      wait_for(syncq);
      for (int i = 0; i < n; ++i) {
         atoms::x[i] = arrx[i];
         atoms::y[i] = arry[i];
         atoms::z[i] = arrz[i];
      }
   }

   {
      std::vector<vel_prec> arrx(n), arry(n), arrz(n);
      darray::copyout(syncq, n, arrx.data(), dup_buf_vx);
      darray::copyout(syncq, n, arry.data(), dup_buf_vy);
      darray::copyout(syncq, n, arrz.data(), dup_buf_vz);
      wait_for(syncq);
      for (int i = 0; i < n; ++i) {
         int j = 3 * i;
         moldyn::v[j] = arrx[i];
         moldyn::v[j + 1] = arry[i];
         moldyn::v[j + 2] = arrz[i];
      }
   }

   {
      std::vector<double> arrx(n), arry(n), arrz(n);
      copy_gradient(calc::grad, arrx.data(), arry.data(), arrz.data(),
                    dup_buf_gx, dup_buf_gy, dup_buf_gz, false);
      // convert gradient to acceleration
      const double ekcal = units::ekcal;
      for (int i = 0; i < n; ++i) {
         int j = 3 * i;
         deriv::desum[j + 0] = arrx[i];
         deriv::desum[j + 1] = arry[i];
         deriv::desum[j + 2] = arrz[i];
         double invmass = 1.0 / atomid::mass[i];
         moldyn::a[j + 0] = -ekcal * arrx[i] * invmass;
         moldyn::a[j + 1] = -ekcal * arry[i] * invmass;
         moldyn::a[j + 2] = -ekcal * arrz[i] * invmass;
         moldyn::aalt[j + 0] = 0;
         moldyn::aalt[j + 1] = 0;
         moldyn::aalt[j + 2] = 0;
      }
   }

   if (mdsave_use_uind()) {
      darray::copyout(syncq, n, polar::uind, dup_buf_uind);
      wait_for(syncq);
   }

   // Record mdsave_end_event when stream (0) is available.
   // nonblk will wait until mdsave_end_event is recorded, so that the dup_
   // arrays are idle and ready to be written.
   dup_stream.end_copyout();

   double dt1 = dt;
   double epot1 = epot;
   double eksum1 = eksum;
   TINKER_RT(mdsave)(&istep, &dt1, &epot1, &eksum1);

   mtx_write.lock();
   idle_write = true;
   cv_write.notify_all();
   mtx_write.unlock();
}
}


void mdsave_async(int istep, time_prec dt)
{
   std::unique_lock<std::mutex> lck_write(mtx_write);
   cv_write.wait(lck_write, [=]() { return idle_write; });
   idle_write = false;


   fut_dup_then_write =
      std::async(std::launch::async, mdsave_dup_then_write, istep, dt);


   std::unique_lock<std::mutex> lck_copy(mtx_dup);
   cv_dup.wait(lck_copy, [=]() { return idle_dup; });
   idle_dup = false;
}


void mdsave_synchronize()
{
   if (fut_dup_then_write.valid())
      fut_dup_then_write.get();
}


void mdsave_data(rc_op op)
{
   if (op & rc_dealloc) {
      dup_stream.deallocate();

      if (mdsave_use_uind()) {
         darray::deallocate(dup_buf_uind);
      }

      darray::deallocate(dup_buf_x, dup_buf_y, dup_buf_z);
      darray::deallocate(dup_buf_vx, dup_buf_vy, dup_buf_vz);
      darray::deallocate(dup_buf_gx, dup_buf_gy, dup_buf_gz);
   }

   if (op & rc_alloc) {
      dup_stream.allocate();

      if (mdsave_use_uind()) {
         darray::allocate(n, &dup_buf_uind);
      } else {
         dup_buf_uind = nullptr;
      }

      darray::allocate(n, &dup_buf_x, &dup_buf_y, &dup_buf_z);
      darray::allocate(n, &dup_buf_vx, &dup_buf_vy, &dup_buf_vz);
      darray::allocate(n, &dup_buf_gx, &dup_buf_gy, &dup_buf_gz);
   }

   if (op & rc_init) {
      idle_dup = false;
      idle_write = true;
   }
}
}
