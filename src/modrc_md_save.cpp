#include "mod_box.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_io.h"
#include "util_potent.h"
#include "util_potential.h"
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
static cudaStream_t dup_stream_uind_;
static real* dup_buf_uind_;
static cudaStream_t dup_stream_bxyz_, dup_stream_v_, dup_stream_g_;
static real* dup_buf_esum_;
static box_t* dup_buf_box_;
static real *dup_buf_x_, *dup_buf_y_, *dup_buf_z_;
static real *dup_buf_vx_, *dup_buf_vy_, *dup_buf_vz_;
static real *dup_buf_gx_, *dup_buf_gy_, *dup_buf_gz_;

void mdsave_data(rc_t rc) {
  if (rc & rc_dealloc) {
    if (mdsave_use_uind_()) {
      check_cudart(cudaStreamDestroy(dup_stream_uind_));
      check_cudart(cudaFree(dup_buf_uind_));
    }

    check_cudart(cudaStreamDestroy(dup_stream_bxyz_));
    check_cudart(cudaStreamDestroy(dup_stream_v_));
    check_cudart(cudaStreamDestroy(dup_stream_g_));

    check_cudart(cudaFree(dup_buf_esum_));
    check_cudart(cudaFree(dup_buf_box_));
    check_cudart(cudaFree(dup_buf_x_));
    check_cudart(cudaFree(dup_buf_y_));
    check_cudart(cudaFree(dup_buf_z_));

    check_cudart(cudaFree(dup_buf_vx_));
    check_cudart(cudaFree(dup_buf_vy_));
    check_cudart(cudaFree(dup_buf_vz_));

    check_cudart(cudaFree(dup_buf_gx_));
    check_cudart(cudaFree(dup_buf_gy_));
    check_cudart(cudaFree(dup_buf_gz_));
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);

    if (mdsave_use_uind_()) {
      check_cudart(cudaStreamCreate(&dup_stream_uind_));
      check_cudart(cudaMalloc(&dup_buf_uind_, rs * 3 * n));
    } else {
      dup_stream_uind_ = nullptr;
      dup_buf_uind_ = nullptr;
    }

    check_cudart(cudaStreamCreate(&dup_stream_bxyz_));
    check_cudart(cudaStreamCreate(&dup_stream_v_));
    check_cudart(cudaStreamCreate(&dup_stream_g_));

    check_cudart(cudaMalloc(&dup_buf_esum_, sizeof(real)));
    check_cudart(cudaMalloc(&dup_buf_box_, sizeof(box_t)));
    check_cudart(cudaMalloc(&dup_buf_x_, rs * n));
    check_cudart(cudaMalloc(&dup_buf_y_, rs * n));
    check_cudart(cudaMalloc(&dup_buf_z_, rs * n));

    check_cudart(cudaMalloc(&dup_buf_vx_, rs * n));
    check_cudart(cudaMalloc(&dup_buf_vy_, rs * n));
    check_cudart(cudaMalloc(&dup_buf_vz_, rs * n));

    check_cudart(cudaMalloc(&dup_buf_gx_, rs * n));
    check_cudart(cudaMalloc(&dup_buf_gy_, rs * n));
    check_cudart(cudaMalloc(&dup_buf_gz_, rs * n));
  }

  if (rc & rc_copyin) {
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
      energy_potential(use_data & calc::vmask);
    }
  }
}

static void mdsave_dup_then_write_(int istep, real dt) {

  // duplicate

  const size_t rs = sizeof(real);

  check_cudart(cudaMemcpyAsync(dup_buf_esum_, esum, rs,
                               cudaMemcpyDeviceToDevice, dup_stream_bxyz_));
  check_cudart(cudaMemcpyAsync(dup_buf_box_, box, sizeof(box_t),
                               cudaMemcpyDeviceToDevice, dup_stream_bxyz_));
  check_cudart(cudaMemcpyAsync(dup_buf_x_, x, rs * n, cudaMemcpyDeviceToDevice,
                               dup_stream_bxyz_));
  check_cudart(cudaMemcpyAsync(dup_buf_y_, y, rs * n, cudaMemcpyDeviceToDevice,
                               dup_stream_bxyz_));
  check_cudart(cudaMemcpyAsync(dup_buf_z_, z, rs * n, cudaMemcpyDeviceToDevice,
                               dup_stream_bxyz_));

  check_cudart(cudaMemcpyAsync(dup_buf_vx_, vx, rs * n,
                               cudaMemcpyDeviceToDevice, dup_stream_v_));
  check_cudart(cudaMemcpyAsync(dup_buf_vy_, vy, rs * n,
                               cudaMemcpyDeviceToDevice, dup_stream_v_));
  check_cudart(cudaMemcpyAsync(dup_buf_vz_, vz, rs * n,
                               cudaMemcpyDeviceToDevice, dup_stream_v_));

  check_cudart(cudaMemcpyAsync(dup_buf_gx_, gx, rs * n,
                               cudaMemcpyDeviceToDevice, dup_stream_g_));
  check_cudart(cudaMemcpyAsync(dup_buf_gy_, gy, rs * n,
                               cudaMemcpyDeviceToDevice, dup_stream_g_));
  check_cudart(cudaMemcpyAsync(dup_buf_gz_, gz, rs * n,
                               cudaMemcpyDeviceToDevice, dup_stream_g_));

  if (mdsave_use_uind_()) {
    check_cudart(cudaMemcpyAsync(dup_buf_uind_, &uind[0][0], rs * 3 * n,
                                 cudaMemcpyDeviceToDevice, dup_stream_uind_));
  }

  check_cudart(cudaStreamSynchronize(dup_stream_bxyz_));
  check_cudart(cudaStreamSynchronize(dup_stream_g_));
  check_cudart(cudaStreamSynchronize(dup_stream_v_));
  if (mdsave_use_uind_())
    check_cudart(cudaStreamSynchronize(dup_stream_uind_));

  mtx_dup.lock();
  idle_dup = true;
  cv_dup.notify_all();
  mtx_dup.unlock();

  // get gpu buffer and write to external files

  box_t b;
  std::vector<real> arrx(n), arry(n), arrz(n);

  check_cudart(cudaMemcpy(&epot, dup_buf_esum_, rs, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(&b, dup_buf_box_, sizeof(box_t), cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arrx.data(), dup_buf_x_, rs * n, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arry.data(), dup_buf_y_, rs * n, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arrz.data(), dup_buf_z_, rs * n, cudaMemcpyDeviceToHost));
  box_data_copyout(b);
  for (int i = 0; i < n; ++i) {
    atoms::x[i] = arrx[i];
    atoms::y[i] = arry[i];
    atoms::z[i] = arrz[i];
  }

  check_cudart(
      cudaMemcpy(arrx.data(), dup_buf_vx_, rs * n, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arry.data(), dup_buf_vy_, rs * n, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arrz.data(), dup_buf_vz_, rs * n, cudaMemcpyDeviceToHost));
  for (int i = 0; i < n; ++i) {
    int j = 3 * i;
    moldyn::v[j] = arrx[i];
    moldyn::v[j + 1] = arry[i];
    moldyn::v[j + 2] = arrz[i];
  }

  check_cudart(
      cudaMemcpy(arrx.data(), dup_buf_gx_, rs * n, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arry.data(), dup_buf_gy_, rs * n, cudaMemcpyDeviceToHost));
  check_cudart(
      cudaMemcpy(arrz.data(), dup_buf_gz_, rs * n, cudaMemcpyDeviceToHost));
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
    check_cudart(cudaMemcpy(arrx.data(), dup_buf_uind_, 3 * rs * n,
                            cudaMemcpyDeviceToHost));
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
