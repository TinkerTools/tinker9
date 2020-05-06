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
   dup_stream.copy_bytes(dup_buf_x, xpos, sizeof(vel_prec) * n);
   dup_stream.copy_bytes(dup_buf_y, ypos, sizeof(vel_prec) * n);
   dup_stream.copy_bytes(dup_buf_z, zpos, sizeof(vel_prec) * n);

   dup_stream.copy_bytes(dup_buf_vx, vx, sizeof(vel_prec) * n);
   dup_stream.copy_bytes(dup_buf_vy, vy, sizeof(vel_prec) * n);
   dup_stream.copy_bytes(dup_buf_vz, vz, sizeof(vel_prec) * n);

   dup_stream.copy_bytes(dup_buf_gx, gx, sizeof(grad_prec) * n);
   dup_stream.copy_bytes(dup_buf_gy, gy, sizeof(grad_prec) * n);
   dup_stream.copy_bytes(dup_buf_gz, gz, sizeof(grad_prec) * n);

   if (mdsave_use_uind()) {
      dup_stream.copy_bytes(&dup_buf_uind[0][0], &uind[0][0],
                            sizeof(real) * 3 * n);
   }

   dup_stream.synchronize();

   mtx_dup.lock();
   idle_dup = true;
   cv_dup.notify_all();
   mtx_dup.unlock();

   // get gpu buffer and write to external files

   energy_prec epot = dup_buf_esum;
   set_tinker_box_module(dup_buf_box);
   if (sizeof(pos_prec) == sizeof(double)) {
      darray::copyout(PROCEED_NEW_Q, n, atoms::x, dup_buf_x);
      darray::copyout(PROCEED_NEW_Q, n, atoms::y, dup_buf_y);
      darray::copyout(WAIT_NEW_Q, n, atoms::z, dup_buf_z);
   } else {
      std::vector<pos_prec> arrx(n), arry(n), arrz(n);
      darray::copyout(PROCEED_NEW_Q, n, arrx.data(), dup_buf_x);
      darray::copyout(PROCEED_NEW_Q, n, arry.data(), dup_buf_y);
      darray::copyout(WAIT_NEW_Q, n, arrz.data(), dup_buf_z);
      for (int i = 0; i < n; ++i) {
         atoms::x[i] = arrx[i];
         atoms::y[i] = arry[i];
         atoms::z[i] = arrz[i];
      }
   }

   {
      std::vector<vel_prec> arrx(n), arry(n), arrz(n);
      darray::copyout(PROCEED_NEW_Q, n, arrx.data(), dup_buf_vx);
      darray::copyout(PROCEED_NEW_Q, n, arry.data(), dup_buf_vy);
      darray::copyout(WAIT_NEW_Q, n, arrz.data(), dup_buf_vz);
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
                    dup_buf_gx, dup_buf_gy, dup_buf_gz);
      // convert gradient to acceleration
      const double ekcal = units::ekcal;
      for (int i = 0; i < n; ++i) {
         int j = 3 * i;
         double invmass = 1.0 / atomid::mass[i];
         moldyn::a[j] = -ekcal * arrx[i] * invmass;
         moldyn::a[j + 1] = -ekcal * arry[i] * invmass;
         moldyn::a[j + 2] = -ekcal * arrz[i] * invmass;
      }
   }

   if (mdsave_use_uind()) {
      darray::copyout(WAIT_NEW_Q, n, polar::uind, dup_buf_uind);
   }

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
