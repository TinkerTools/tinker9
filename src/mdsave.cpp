#include "ff/amoebamod.h"
#include "ff/energy.h"
#include "ff/potent.h"
#include "md/misc.h"
#include "md/pq.h"
#include "tool/execq.h"
#include "tool/iofortstr.h"
#include <condition_variable>
#include <future>
#include <mutex>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/deriv.hh>
#include <tinker/detail/files.hh>
#include <tinker/detail/moldyn.hh>
#include <tinker/detail/output.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/titles.hh>
#include <tinker/detail/units.hh>
#include <tinker/routines.h>

#if TINKER_CUDART
#   include "tool/error.h"
#   include "tool/gpucard.h"
#   include <cuda_runtime.h>
#endif

namespace tinker {
static std::mutex mtx_dup, mtx_write;
static std::condition_variable cv_dup, cv_write;
static bool idle_dup, idle_write;
static std::future<void> fut_dup_then_write;

static bool mdsaveUseUind()
{
   return static_cast<bool>(output::uindsave) and use(Potent::POLAR);
}

static ExecQ dup_stream;
static real (*dup_buf_uind)[3];
static energy_prec dup_buf_esum;
static Box dup_buf_box;
static pos_prec *dup_buf_x, *dup_buf_y, *dup_buf_z;
static vel_prec *dup_buf_vx, *dup_buf_vy, *dup_buf_vz;
static grad_prec *dup_buf_gx, *dup_buf_gy, *dup_buf_gz;

static void mdsaveDupThenWrite(int istep, time_prec dt)
{
#if TINKER_CUDART
   // This function (mdsaveDupThenWrite) will run in another CPU thread.
   // There is no guarantee that the CUDA runtime will use the same GPU card as
   // the main thread, unless cudaSetDevice() is called explicitly.
   //
   // Of course this is not a problem if the computer has only one GPU card.
   check_rt(cudaSetDevice(idevice));
#endif

   // duplicate

   dup_buf_esum = esum;
   boxGetCurrent(dup_buf_box);
   darray::copy(g::q0, n, dup_buf_x, xpos);
   darray::copy(g::q0, n, dup_buf_y, ypos);
   darray::copy(g::q0, n, dup_buf_z, zpos);

   darray::copy(g::q0, n, dup_buf_vx, vx);
   darray::copy(g::q0, n, dup_buf_vy, vy);
   darray::copy(g::q0, n, dup_buf_vz, vz);

   darray::copy(g::q0, n, dup_buf_gx, gx);
   darray::copy(g::q0, n, dup_buf_gy, gy);
   darray::copy(g::q0, n, dup_buf_gz, gz);

   if (mdsaveUseUind()) {
      darray::copy(g::q0, 3 * n, &dup_buf_uind[0][0], &uind[0][0]);
   }

   // Record mdsave_begin_event when g::s0 is available.
   // g::s1 will wait until mdsave_begin_event is recorded.
   dup_stream.beginCopyout();

   mtx_dup.lock();
   idle_dup = true;
   cv_dup.notify_all();
   mtx_dup.unlock();

   // get gpu buffer and write to external files

   energy_prec epot = dup_buf_esum;
   boxSetTinker(dup_buf_box);
   if (sizeof(pos_prec) == sizeof(double)) {
      darray::copyout(g::q1, n, atoms::x, dup_buf_x);
      darray::copyout(g::q1, n, atoms::y, dup_buf_y);
      darray::copyout(g::q1, n, atoms::z, dup_buf_z);
      waitFor(g::q1);
   } else {
      std::vector<pos_prec> arrx(n), arry(n), arrz(n);
      darray::copyout(g::q1, n, arrx.data(), dup_buf_x);
      darray::copyout(g::q1, n, arry.data(), dup_buf_y);
      darray::copyout(g::q1, n, arrz.data(), dup_buf_z);
      waitFor(g::q1);
      for (int i = 0; i < n; ++i) {
         atoms::x[i] = arrx[i];
         atoms::y[i] = arry[i];
         atoms::z[i] = arrz[i];
      }
   }

   {
      std::vector<vel_prec> arrx(n), arry(n), arrz(n);
      darray::copyout(g::q1, n, arrx.data(), dup_buf_vx);
      darray::copyout(g::q1, n, arry.data(), dup_buf_vy);
      darray::copyout(g::q1, n, arrz.data(), dup_buf_vz);
      waitFor(g::q1);
      for (int i = 0; i < n; ++i) {
         int j = 3 * i;
         moldyn::v[j] = arrx[i];
         moldyn::v[j + 1] = arry[i];
         moldyn::v[j + 2] = arrz[i];
      }
   }

   {
      std::vector<double> arrx(n), arry(n), arrz(n);
      copyGradientSync(calc::grad, arrx.data(), arry.data(), arrz.data(), dup_buf_gx, dup_buf_gy,
         dup_buf_gz, g::q1);
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

   if (mdsaveUseUind()) {
      darray::copyout(g::q1, n, polar::uind, dup_buf_uind);
      waitFor(g::q1);
   }

   // Record mdsave_end_event when g::s1 is available.
   // g::s0 will wait until mdsave_end_event is recorded, so that the dup_
   // arrays are idle and ready to be written.
   dup_stream.endCopyout();

   double dt1 = dt;
   double epot1 = epot;
   double eksum1 = eksum;
   tinker_f_mdsave(&istep, &dt1, &epot1, &eksum1);

   mtx_write.lock();
   idle_write = true;
   cv_write.notify_all();
   mtx_write.unlock();
}
}

namespace tinker {
void mdsaveAsync(int istep, time_prec dt)
{
   std::unique_lock<std::mutex> lck_write(mtx_write);
   cv_write.wait(lck_write, [=]() { return idle_write; });
   idle_write = false;

   fut_dup_then_write = std::async(std::launch::async, mdsaveDupThenWrite, istep, dt);

   std::unique_lock<std::mutex> lck_copy(mtx_dup);
   cv_dup.wait(lck_copy, [=]() { return idle_dup; });
   idle_dup = false;
}

void mdDebugSaveSync()
{
   auto DOT3 = [](const double* a, const double* b) -> double {
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   };

   static bool first = true;
   static std::vector<pos_prec> qx, qy, qz;
   static std::string fname, title;
   if (first) {
      first = false;
      qx.resize(n);
      qy.resize(n);
      qz.resize(n);
      FstrView fstr = files::filename;
      fname = fstr.trim();
      auto pos = fname.find_last_of('.');
      fname = fname.substr(0, pos) + ".dbg";
      FstrView ftitl = titles::title;
      title = ftitl.trim();
      FILE* ftmp = fopen(fname.c_str(), "w");
      fclose(ftmp);
   }

   FILE* fout = fopen(fname.c_str(), "a");

   darray::copyout(g::q0, n, qx.data(), xpos);
   darray::copyout(g::q0, n, qy.data(), ypos);
   darray::copyout(g::q0, n, qz.data(), zpos);
   waitFor(g::q0);

   bool bign = n > 999999;
   if (bign)
      fprintf(fout, "%8d  %s\n", n, title.c_str());
   else
      fprintf(fout, "%6d  %s\n", n, title.c_str());

   if (box_shape != BoxShape::UNBOUND) {
      double ax[3] = {lvec1.x, lvec2.x, lvec3.x};
      double bx[3] = {lvec1.y, lvec2.y, lvec3.y};
      double cx[3] = {lvec1.z, lvec2.z, lvec3.z};

      double xb = std::sqrt(DOT3(ax, ax));
      double yb = std::sqrt(DOT3(bx, bx));
      double zb = std::sqrt(DOT3(cx, cx));

      double cos_a = DOT3(bx, cx) / (yb * zb);
      double cos_b = DOT3(cx, ax) / (zb * xb);
      double cos_c = DOT3(ax, bx) / (xb * yb);

      double al = 90.0;
      double be = 90.0;
      double ga = 90.0;
      if (cos_a != 0.0) al = (180 / M_PI) * std::acos(cos_a);
      if (cos_b != 0.0) be = (180 / M_PI) * std::acos(cos_b);
      if (cos_c != 0.0) ga = (180 / M_PI) * std::acos(cos_c);

      if (bign)
         fprintf(fout, " %14.6lf%14.6lf%14.6lf%14.6lf%14.6lf%14.6lf\n", xb, yb, zb, al, be, ga);
      else
         fprintf(fout, " %12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf\n", xb, yb, zb, al, be, ga);
   }

   const char* fcord;
   const char* fcoup;
   if (bign) {
      fcord = "%8d  %c%c%c%14.6lf%14.6lf%14.6lf%8d";
      fcoup = "%8d";
   } else {
      fcord = "%6d  %c%c%c%12.6lf%12.6lf%12.6lf%6d";
      fcoup = "%6d";
   }
   for (int i = 0; i < n; ++i) {
      const char* nm = atomid::name[i];
      fprintf(fout, fcord, i + 1, nm[0], nm[1], nm[2], qx[i], qy[i], qz[i], atoms::type[i]);
      for (int k = 0; k < couple::n12[i]; ++k)
         fprintf(fout, fcoup, couple::i12[i][k]);
      fprintf(fout, "%s", "\n");
   }

   fclose(fout);
}

void mdsaveSynchronize()
{
   if (fut_dup_then_write.valid())
      fut_dup_then_write.get();
}

void mdsaveData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
      dup_stream.deallocate();

      if (mdsaveUseUind()) {
         darray::deallocate(dup_buf_uind);
      }

      darray::deallocate(dup_buf_x, dup_buf_y, dup_buf_z);
      darray::deallocate(dup_buf_vx, dup_buf_vy, dup_buf_vz);
      darray::deallocate(dup_buf_gx, dup_buf_gy, dup_buf_gz);
   }

   if (op & RcOp::ALLOC) {
      dup_stream.allocate();

      if (mdsaveUseUind()) {
         darray::allocate(n, &dup_buf_uind);
      } else {
         dup_buf_uind = nullptr;
      }

      darray::allocate(n, &dup_buf_x, &dup_buf_y, &dup_buf_z);
      darray::allocate(n, &dup_buf_vx, &dup_buf_vy, &dup_buf_vz);
      darray::allocate(n, &dup_buf_gx, &dup_buf_gy, &dup_buf_gz);
   }

   if (op & RcOp::INIT) {
      idle_dup = false;
      idle_write = true;
   }
}
}
