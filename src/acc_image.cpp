#include "acc_seq.h"

TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
static inline void box_null_image() {}

#pragma acc routine seq
static inline void box_null_imagen() {}

#pragma acc routine seq
static inline void box_ortho_image(real& __restrict__ xr, real& __restrict__ yr,
                                   real& __restrict__ zr,
                                   const Box* __restrict__ pb) {
  real fx = xr * pb->recip[0][0];
  real fy = yr * pb->recip[1][1];
  real fz = zr * pb->recip[2][2];
  fx = REAL_FLOOR(fx + 0.5f);
  fy = REAL_FLOOR(fy + 0.5f);
  fz = REAL_FLOOR(fz + 0.5f);
  xr -= (fx * pb->lvec[0][0]);
  yr -= (fy * pb->lvec[1][1]);
  zr -= (fz * pb->lvec[2][2]);
}

#pragma acc routine seq
static inline void box_ortho_imagen(real& __restrict__ xr,
                                    real& __restrict__ yr,
                                    real& __restrict__ zr,
                                    const Box* __restrict__ pb) {
  real lx = pb->lvec[0][0];
  real lx2 = lx * 0.5f;
  real a = REAL_ABS(xr);
  xr = (a > lx2 ? a - lx : a);

  real ly = pb->lvec[1][1];
  real ly2 = ly * 0.5f;
  real b = REAL_ABS(yr);
  yr = (b > ly2 ? b - ly : b);

  real lz = pb->lvec[2][2];
  real lz2 = lz * 0.5f;
  real c = REAL_ABS(zr);
  zr = (c > lz2 ? c - lz : c);
}

#pragma acc routine seq
static inline void box_mono_image(real& __restrict__ xr, real& __restrict__ yr,
                                  real& __restrict__ zr,
                                  const Box* __restrict__ pb) {
  real fx = xr * pb->recip[0][0] + zr * pb->recip[0][2];
  real fy = yr * pb->recip[1][1];
  real fz = zr * pb->recip[2][2];
  fx = REAL_FLOOR(fx + 0.5f);
  fy = REAL_FLOOR(fy + 0.5f);
  fz = REAL_FLOOR(fz + 0.5f);
  xr -= (fx * pb->lvec[0][0] + fz * pb->lvec[0][2]);
  yr -= (fy * pb->lvec[1][1]);
  zr -= (fz * pb->lvec[2][2]);
}

#pragma acc routine seq
static inline void box_mono_imagen(real& __restrict__ xr, real& __restrict__ yr,
                                   real& __restrict__ zr,
                                   const Box* __restrict__ pb) {
  // TODO: a real imagen routine
  box_mono_image(xr, yr, zr, pb);
}

#pragma acc routine seq
static inline void box_tri_image(real& __restrict__ xr, real& __restrict__ yr,
                                 real& __restrict__ zr,
                                 const Box* __restrict__ pb) {
  real fx = xr * pb->recip[0][0] + yr * pb->recip[0][1] + zr * pb->recip[0][2];
  real fy = yr * pb->recip[1][1] + zr * pb->recip[1][2];
  real fz = zr * pb->recip[2][2];
  fx = REAL_FLOOR(fx + 0.5f);
  fy = REAL_FLOOR(fy + 0.5f);
  fz = REAL_FLOOR(fz + 0.5f);
  xr -= (fx * pb->lvec[0][0] + fy * pb->lvec[0][1] + fz * pb->lvec[0][2]);
  yr -= (fy * pb->lvec[1][1] + fz * pb->lvec[1][2]);
  zr -= (fz * pb->lvec[2][2]);
}

#pragma acc routine seq
static inline void box_tri_imagen(real& __restrict__ xr, real& __restrict__ yr,
                                  real& __restrict__ zr,
                                  const Box* __restrict__ pb) {
  // TODO: a real imagen routine
  box_tri_image(xr, yr, zr, pb);
}

#pragma acc routine seq
void image(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
           const Box* __restrict__ pb) {
  switch (pb->shape) {
  case Box::ortho:
    box_ortho_image(xr, yr, zr, pb);
    break;
  case Box::mono:
    box_mono_image(xr, yr, zr, pb);
    break;
  case Box::tri:
    box_tri_image(xr, yr, zr, pb);
    break;
  default:
    // Box::null
    box_null_image();
    break;
  }
}

#pragma acc routine seq
void imagen(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
            const Box* __restrict__ pb) {
  switch (pb->shape) {
  case Box::ortho:
    box_ortho_imagen(xr, yr, zr, pb);
    break;
  case Box::mono:
    box_mono_imagen(xr, yr, zr, pb);
    break;
  case Box::tri:
    box_tri_imagen(xr, yr, zr, pb);
    break;
  default:
    // Box::null
    box_null_imagen();
    break;
  }
}
TINKER_NAMESPACE_END
