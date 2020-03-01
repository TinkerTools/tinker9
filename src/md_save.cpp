#include "box.h"
#include "e_polar.h"
#include "energy.h"
#include "execq.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include "tinker_rt.h"
#include <condition_variable>
#include <fstream>
#include <future>
#include <mutex>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/files.hh>
#include <tinker/detail/moldyn.hh>
#include <tinker/detail/output.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/units.hh>

TINKER_NAMESPACE_BEGIN
static std::mutex mtx_dup, mtx_write;
static std::condition_variable cv_dup, cv_write;
static bool idle_dup, idle_write;
static std::future<void> fut_dup_then_write;

static bool mdsave_use_uind_()
{
   return output::uindsave && use_potent(polar_term);
}
static ExecQ dup_stream_;
static real (*dup_buf_uind_)[3];
static energy_prec dup_buf_esum_;
static Box dup_buf_box_;
static pos_prec *dup_buf_x_, *dup_buf_y_, *dup_buf_z_;
static vel_prec *dup_buf_vx_, *dup_buf_vy_, *dup_buf_vz_;
static grad_prec *dup_buf_gx_, *dup_buf_gy_, *dup_buf_gz_;

void mdsave_data(rc_op op)
{
   if (op & rc_dealloc) {
      dup_stream_.deallocate();

      if (mdsave_use_uind_()) {
         device_array::deallocate(dup_buf_uind_);
      }

      device_array::deallocate(dup_buf_x_, dup_buf_y_, dup_buf_z_);
      device_array::deallocate(dup_buf_vx_, dup_buf_vy_, dup_buf_vz_);
      device_array::deallocate(dup_buf_gx_, dup_buf_gy_, dup_buf_gz_);
   }

   if (op & rc_alloc) {
      dup_stream_.allocate();

      if (mdsave_use_uind_()) {
         device_array::allocate(n, &dup_buf_uind_);
      } else {
         dup_buf_uind_ = nullptr;
      }

      device_array::allocate(n, &dup_buf_x_, &dup_buf_y_, &dup_buf_z_);
      device_array::allocate(n, &dup_buf_vx_, &dup_buf_vy_, &dup_buf_vz_);
      device_array::allocate(n, &dup_buf_gx_, &dup_buf_gy_, &dup_buf_gz_);
   }

   if (op & rc_init) {
      idle_dup = false;
      idle_write = true;
   }
}

static void mdsave_dup_then_write_(int istep, time_prec dt)
{

   // duplicate

   dup_buf_esum_ = esum;
   get_default_box(dup_buf_box_);
   dup_stream_.copy_bytes(dup_buf_x_, xpos, sizeof(vel_prec) * n);
   dup_stream_.copy_bytes(dup_buf_y_, ypos, sizeof(vel_prec) * n);
   dup_stream_.copy_bytes(dup_buf_z_, zpos, sizeof(vel_prec) * n);

   dup_stream_.copy_bytes(dup_buf_vx_, vx, sizeof(vel_prec) * n);
   dup_stream_.copy_bytes(dup_buf_vy_, vy, sizeof(vel_prec) * n);
   dup_stream_.copy_bytes(dup_buf_vz_, vz, sizeof(vel_prec) * n);

   dup_stream_.copy_bytes(dup_buf_gx_, gx, sizeof(grad_prec) * n);
   dup_stream_.copy_bytes(dup_buf_gy_, gy, sizeof(grad_prec) * n);
   dup_stream_.copy_bytes(dup_buf_gz_, gz, sizeof(grad_prec) * n);

   if (mdsave_use_uind_()) {
      dup_stream_.copy_bytes(&dup_buf_uind_[0][0], &uind[0][0],
                             sizeof(real) * 3 * n);
   }

   dup_stream_.synchronize();

   mtx_dup.lock();
   idle_dup = true;
   cv_dup.notify_all();
   mtx_dup.unlock();

   // get gpu buffer and write to external files


   energy_prec epot = dup_buf_esum_;
   set_tinker_box_module(dup_buf_box_);
   if (sizeof(pos_prec) == sizeof(double)) {
      device_array::copyout(PROCEED_NEW_Q, n, atoms::x, dup_buf_x_);
      device_array::copyout(PROCEED_NEW_Q, n, atoms::y, dup_buf_y_);
      device_array::copyout(WAIT_NEW_Q, n, atoms::z, dup_buf_z_);
   } else {
      std::vector<pos_prec> arrx(n), arry(n), arrz(n);
      device_array::copyout(PROCEED_NEW_Q, n, arrx.data(), dup_buf_x_);
      device_array::copyout(PROCEED_NEW_Q, n, arry.data(), dup_buf_y_);
      device_array::copyout(WAIT_NEW_Q, n, arrz.data(), dup_buf_z_);
      for (int i = 0; i < n; ++i) {
         atoms::x[i] = arrx[i];
         atoms::y[i] = arry[i];
         atoms::z[i] = arrz[i];
      }
   }

   {
      std::vector<vel_prec> arrx(n), arry(n), arrz(n);
      device_array::copyout(PROCEED_NEW_Q, n, arrx.data(), dup_buf_vx_);
      device_array::copyout(PROCEED_NEW_Q, n, arry.data(), dup_buf_vy_);
      device_array::copyout(WAIT_NEW_Q, n, arrz.data(), dup_buf_vz_);
      for (int i = 0; i < n; ++i) {
         int j = 3 * i;
         moldyn::v[j] = arrx[i];
         moldyn::v[j + 1] = arry[i];
         moldyn::v[j + 2] = arrz[i];
      }
   }

   {
      std::vector<double> arrx(n), arry(n), arrz(n);
      copy_energy(calc::grad, nullptr, arrx.data(), arry.data(), arrz.data(),
                  nullptr, dup_buf_gx_, dup_buf_gy_, dup_buf_gz_);
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

   if (mdsave_use_uind_()) {
      device_array::copyout(WAIT_NEW_Q, n, polar::uind, dup_buf_uind_);
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

void mdsave_async(int istep, time_prec dt)
{
   std::unique_lock<std::mutex> lck_write(mtx_write);
   cv_write.wait(lck_write, [=]() { return idle_write; });
   idle_write = false;

   fut_dup_then_write =
      std::async(std::launch::async, mdsave_dup_then_write_, istep, dt);

   std::unique_lock<std::mutex> lck_copy(mtx_dup);
   cv_dup.wait(lck_copy, [=]() { return idle_dup; });
   idle_dup = false;
}

void mdsave_synchronize()
{
   if (fut_dup_then_write.valid())
      fut_dup_then_write.get();
}
TINKER_NAMESPACE_END
