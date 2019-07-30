#include "util_elec.h"
#include "util_mdstate.h"

TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
static void rotsite(int isite, const real (*__restrict__ a)[3],
                    real (*__restrict__ rpole)[10],
                    const real (*__restrict__ pole)[10]) {
  static_assert(mpl_total == 10, "");

  // charge
  rpole[isite][0] = pole[isite][0];
  // dipole
  for (int i = 1; i < 4; ++i) {
    rpole[isite][i] = 0;
    for (int j = 1; j < 4; ++j)
      rpole[isite][i] += pole[isite][j] * a[j - 1][i - 1];
  }
  // quadrupole
  real rp[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  real mp[3][3];
  mp[_x][_x] = pole[isite][mpl_pme_xx];
  mp[_x][_y] = pole[isite][mpl_pme_xy];
  mp[_x][_z] = pole[isite][mpl_pme_xz];
  mp[_y][_x] = pole[isite][mpl_pme_yx];
  mp[_y][_y] = pole[isite][mpl_pme_yy];
  mp[_y][_z] = pole[isite][mpl_pme_yz];
  mp[_z][_x] = pole[isite][mpl_pme_zx];
  mp[_z][_y] = pole[isite][mpl_pme_zy];
  mp[_z][_z] = pole[isite][mpl_pme_zz];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // if (j < i) {
      //  rp[j][i] = rp[i][j];
      // } else {
      for (int k = 0; k < 3; ++k)
        for (int m = 0; m < 3; ++m)
          rp[j][i] += a[k][i] * a[m][j] * mp[m][k];
      // }
    }
  }
  rpole[isite][mpl_pme_xx] = rp[_x][_x];
  rpole[isite][mpl_pme_xy] = rp[_x][_y];
  rpole[isite][mpl_pme_xz] = rp[_x][_z];
  rpole[isite][mpl_pme_yy] = rp[_y][_y];
  rpole[isite][mpl_pme_yz] = rp[_y][_z];
  rpole[isite][mpl_pme_zz] = rp[_z][_z];
}

void rotpole() {
  #pragma acc data deviceptr(x,y,z,zaxis,rpole,pole)
  #pragma acc parallel loop
  for (int i = 0; i < n; ++i) {
    // rotmat routine
    real xi = x[i];
    real yi = y[i];
    real zi = z[i];
    int iz = zaxis[i].zaxis;
    int ix = zaxis[i].xaxis;
    int iy = INT_ABS(zaxis[i].yaxis) - 1;
    int polaxe = zaxis[i].polaxe;
    real a[3][3] = {{1, 0, 0}, {0, 0, 0}, {0, 0, 1}};

    if (polaxe == pole_z_only) {
      real dx = x[iz] - xi;
      real dy = y[iz] - yi;
      real dz = z[iz] - zi;
      real r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[2][0] = dx * r1;
      a[2][1] = dy * r1;
      a[2][2] = dz * r1;

      dx = 1;
      dy = 0;
      dz = 0;
      real dot = a[2][0];
      if (REAL_ABS(dot) > (real)0.866) {
        dx = 0;
        dy = 1;
        dot = a[2][1];
      }
      dx -= dot * a[2][0];
      dy -= dot * a[2][1];
      dz -= dot * a[2][2];
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[0][0] = dx * r1;
      a[0][1] = dy * r1;
      a[0][2] = dz * r1;
    } else if (polaxe == pole_z_then_x) {
      real dx = x[iz] - xi;
      real dy = y[iz] - yi;
      real dz = z[iz] - zi;
      real r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[2][0] = dx * r1;
      a[2][1] = dy * r1;
      a[2][2] = dz * r1;

      dx = x[ix] - xi;
      dy = y[ix] - yi;
      dz = z[ix] - zi;
      real dot = dx * a[2][0] + dy * a[2][1] + dz * a[2][2];
      dx -= dot * a[2][0];
      dy -= dot * a[2][1];
      dz -= dot * a[2][2];
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[0][0] = dx * r1;
      a[0][1] = dy * r1;
      a[0][2] = dz * r1;
    } else if (polaxe == pole_bisector) {
      real dx = x[iz] - xi;
      real dy = y[iz] - yi;
      real dz = z[iz] - zi;
      real r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx1 = dx * r1;
      real dy1 = dy * r1;
      real dz1 = dz * r1;

      dx = x[ix] - xi;
      dy = y[ix] - yi;
      dz = z[ix] - zi;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx2 = dx * r1;
      real dy2 = dy * r1;
      real dz2 = dz * r1;

      dx = dx1 + dx2;
      dy = dy1 + dy2;
      dz = dz1 + dz2;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[2][0] = dx * r1;
      a[2][1] = dy * r1;
      a[2][2] = dz * r1;

      real dot = dx2 * a[2][0] + dy2 * a[2][1] + dz2 * a[2][2];
      dx = dx2 - dot * a[2][0];
      dy = dy2 - dot * a[2][1];
      dz = dz2 - dot * a[2][2];
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[0][0] = dx * r1;
      a[0][1] = dy * r1;
      a[0][2] = dz * r1;
    } else if (polaxe == pole_z_bisect) {
      real dx = x[iz] - xi;
      real dy = y[iz] - yi;
      real dz = z[iz] - zi;
      real r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[2][0] = dx * r1;
      a[2][1] = dy * r1;
      a[2][2] = dz * r1;

      dx = x[ix] - xi;
      dy = y[ix] - yi;
      dz = z[ix] - zi;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx1 = dx * r1;
      real dy1 = dy * r1;
      real dz1 = dz * r1;
      dx = x[iy] - xi;
      dy = y[iy] - yi;
      dz = z[iy] - zi;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx2 = dx * r1;
      real dy2 = dy * r1;
      real dz2 = dz * r1;
      dx = dx1 + dx2;
      dy = dy1 + dy2;
      dz = dz1 + dz2;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      dx = dx * r1;
      dy = dy * r1;
      dz = dz * r1;
      real dot = dx * a[2][0] + dy * a[2][1] + dz * a[2][2];
      dx -= dot * a[2][0];
      dy -= dot * a[2][1];
      dz -= dot * a[2][2];
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[0][0] = dx * r1;
      a[0][1] = dy * r1;
      a[0][2] = dz * r1;
    } else if (polaxe == pole_3_fold) {
      real dx = x[iz] - xi;
      real dy = y[iz] - yi;
      real dz = z[iz] - zi;
      real r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx1 = dx * r1;
      real dy1 = dy * r1;
      real dz1 = dz * r1;

      dx = x[ix] - xi;
      dy = y[ix] - yi;
      dz = z[ix] - zi;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx2 = dx * r1;
      real dy2 = dy * r1;
      real dz2 = dz * r1;

      dx = x[iy] - xi;
      dy = y[iy] - yi;
      dz = z[iy] - zi;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      real dx3 = dx * r1;
      real dy3 = dy * r1;
      real dz3 = dz * r1;

      dx = dx1 + dx2 + dx3;
      dy = dy1 + dy2 + dy3;
      dz = dz1 + dz2 + dz3;
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[2][0] = dx * r1;
      a[2][1] = dy * r1;
      a[2][2] = dz * r1;

      real dot = dx2 * a[2][0] + dy2 * a[2][1] + dz2 * a[2][2];
      dx = dx2 - dot * a[2][0];
      dy = dy2 - dot * a[2][1];
      dz = dz2 - dot * a[2][2];
      r1 = REAL_RSQRT(dx * dx + dy * dy + dz * dz);
      a[0][0] = dx * r1;
      a[0][1] = dy * r1;
      a[0][2] = dz * r1;
    }
    a[1][0] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
    a[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
    a[1][2] = a[0][1] * a[2][0] - a[0][0] * a[2][1];

    // rotsite routine
    rotsite(i, a, rpole, pole);
  }
}
TINKER_NAMESPACE_END
