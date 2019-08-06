#include "array.h"
#include "box.h"
#include "energy.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

#include <condition_variable>
#include <fstream>
#include <future>
#include <mutex>

TINKER_NAMESPACE_BEGIN
static std::mutex mtx_dup, mtx_write;
static std::condition_variable cv_dup, cv_write;
static bool idle_dup, idle_write;
static std::future<void> fut_dup_then_write;

static bool mdsave_use_uind_() {
  return output::uindsave && use_potent(polar_term);
}
static Stream dup_stream_uind_;
static real* dup_buf_uind_;
static Stream dup_stream_bxyz_, dup_stream_v_, dup_stream_g_;
static real* dup_buf_esum_;
static Box* dup_buf_box_;
static real *dup_buf_x_, *dup_buf_y_, *dup_buf_z_;
static real *dup_buf_vx_, *dup_buf_vy_, *dup_buf_vz_;
static real *dup_buf_gx_, *dup_buf_gy_, *dup_buf_gz_;

void mdsave_data(rc_op op) {
  if (op & rc_dealloc) {
    if (mdsave_use_uind_()) {
      dealloc_stream(dup_stream_uind_);
      dealloc_bytes(dup_buf_uind_);
    }

    dealloc_stream(dup_stream_bxyz_);
    dealloc_stream(dup_stream_v_);
    dealloc_stream(dup_stream_g_);

    dealloc_bytes(dup_buf_esum_);
    dealloc_bytes(dup_buf_box_);
    dealloc_bytes(dup_buf_x_);
    dealloc_bytes(dup_buf_y_);
    dealloc_bytes(dup_buf_z_);

    dealloc_bytes(dup_buf_vx_);
    dealloc_bytes(dup_buf_vy_);
    dealloc_bytes(dup_buf_vz_);

    dealloc_bytes(dup_buf_gx_);
    dealloc_bytes(dup_buf_gy_);
    dealloc_bytes(dup_buf_gz_);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);

    if (mdsave_use_uind_()) {
      alloc_stream(&dup_stream_uind_);
      alloc_bytes(&dup_buf_uind_, rs * 3 * n);
    } else {
      dup_stream_uind_ = nullptr;
      dup_buf_uind_ = nullptr;
    }

    alloc_stream(&dup_stream_bxyz_);
    alloc_stream(&dup_stream_v_);
    alloc_stream(&dup_stream_g_);

    alloc_bytes(&dup_buf_esum_, sizeof(real));
    alloc_bytes(&dup_buf_box_, sizeof(Box));
    alloc_bytes(&dup_buf_x_, rs * n);
    alloc_bytes(&dup_buf_y_, rs * n);
    alloc_bytes(&dup_buf_z_, rs * n);

    alloc_bytes(&dup_buf_vx_, rs * n);
    alloc_bytes(&dup_buf_vy_, rs * n);
    alloc_bytes(&dup_buf_vz_, rs * n);

    alloc_bytes(&dup_buf_gx_, rs * n);
    alloc_bytes(&dup_buf_gy_, rs * n);
    alloc_bytes(&dup_buf_gz_, rs * n);
  }

  if (op & rc_init) {
    idle_dup = false;
    idle_write = true;

    fstr_view f1 = files::filename;
    fstr_view f2 = f1(1, files::leng);
    std::string dynfile = f2.trim() + ".dyn";
    if (std::ifstream(dynfile)) {

      // convert acceleration to gradient

      std::vector<double> gbuf(n);
      for (int i = 0; i < n; ++i)
        gbuf[i] = -moldyn::a[3 * i] * atomid::mass[i] / units::ekcal;
      copyin_array(gx, gbuf.data(), n);
      for (int i = 0; i < n; ++i)
        gbuf[i] = -moldyn::a[3 * i + 1] * atomid::mass[i] / units::ekcal;
      copyin_array(gy, gbuf.data(), n);
      for (int i = 0; i < n; ++i)
        gbuf[i] = -moldyn::a[3 * i + 2] * atomid::mass[i] / units::ekcal;
      copyin_array(gz, gbuf.data(), n);
    } else {
      energy_potential(rc_flag & calc::vmask);
    }
  }
}

static void mdsave_dup_then_write_(int istep, real dt) {

  // duplicate

  const size_t rs = sizeof(real);

  copy_bytes_async(dup_buf_esum_, esum, rs, dup_stream_bxyz_);
  copy_bytes_async(dup_buf_box_, box, sizeof(Box), dup_stream_bxyz_);
  copy_bytes_async(dup_buf_x_, x, rs * n, dup_stream_bxyz_);
  copy_bytes_async(dup_buf_y_, y, rs * n, dup_stream_bxyz_);
  copy_bytes_async(dup_buf_z_, z, rs * n, dup_stream_bxyz_);

  copy_bytes_async(dup_buf_vx_, vx, rs * n, dup_stream_v_);
  copy_bytes_async(dup_buf_vy_, vy, rs * n, dup_stream_v_);
  copy_bytes_async(dup_buf_vz_, vz, rs * n, dup_stream_v_);

  copy_bytes_async(dup_buf_gx_, gx, rs * n, dup_stream_g_);
  copy_bytes_async(dup_buf_gy_, gy, rs * n, dup_stream_g_);
  copy_bytes_async(dup_buf_gz_, gz, rs * n, dup_stream_g_);

  if (mdsave_use_uind_()) {
    copy_bytes_async(dup_buf_uind_, &uind[0][0], rs * 3 * n, dup_stream_uind_);
  }

  sync_stream(dup_stream_bxyz_);
  sync_stream(dup_stream_g_);
  sync_stream(dup_stream_v_);
  if (mdsave_use_uind_())
    sync_stream(dup_stream_uind_);

  mtx_dup.lock();
  idle_dup = true;
  cv_dup.notify_all();
  mtx_dup.unlock();

  // get gpu buffer and write to external files

  std::vector<real> arrx(n), arry(n), arrz(n);

  copyout_bytes(&epot, dup_buf_esum_, rs);
  copyout_box_data(dup_buf_box_);
  copyout_bytes(arrx.data(), dup_buf_x_, rs * n);
  copyout_bytes(arry.data(), dup_buf_y_, rs * n);
  copyout_bytes(arrz.data(), dup_buf_z_, rs * n);
  for (int i = 0; i < n; ++i) {
    atoms::x[i] = arrx[i];
    atoms::y[i] = arry[i];
    atoms::z[i] = arrz[i];
  }

  copyout_bytes(arrx.data(), dup_buf_vx_, rs * n);
  copyout_bytes(arry.data(), dup_buf_vy_, rs * n);
  copyout_bytes(arrz.data(), dup_buf_vz_, rs * n);
  for (int i = 0; i < n; ++i) {
    int j = 3 * i;
    moldyn::v[j] = arrx[i];
    moldyn::v[j + 1] = arry[i];
    moldyn::v[j + 2] = arrz[i];
  }

  copyout_bytes(arrx.data(), dup_buf_gx_, rs * n);
  copyout_bytes(arry.data(), dup_buf_gy_, rs * n);
  copyout_bytes(arrz.data(), dup_buf_gz_, rs * n);
  // convert gradient to acceleration
  const double ekcal = units::ekcal;
  for (int i = 0; i < n; ++i) {
    int j = 3 * i;
    double invmass = 1.0 / atomid::mass[i];
    moldyn::a[j] = -ekcal * arrx[i] * invmass;
    moldyn::a[j + 1] = -ekcal * arry[i] * invmass;
    moldyn::a[j + 2] = -ekcal * arrz[i] * invmass;
  }

  if (mdsave_use_uind_()) {
    arrx.resize(3 * n);
    copyout_bytes(arrx.data(), dup_buf_uind_, 3 * rs * n);
    for (int i = 0; i < 3 * n; ++i)
      polar::uind[i] = arrx[i];
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

void mdsave_async(int istep, real dt) {
  std::unique_lock<std::mutex> lck_write(mtx_write);
  cv_write.wait(lck_write, [=]() { return idle_write; });
  idle_write = false;

  fut_dup_then_write =
      std::async(std::launch::async, mdsave_dup_then_write_, istep, dt);

  std::unique_lock<std::mutex> lck_copy(mtx_dup);
  cv_dup.wait(lck_copy, [=]() { return idle_dup; });
  idle_dup = false;
}

void mdsave_synchronize() {
  if (fut_dup_then_write.valid())
    fut_dup_then_write.get();
}
TINKER_NAMESPACE_END
