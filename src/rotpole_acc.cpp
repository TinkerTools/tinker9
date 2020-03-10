#include "empole.h"
#include "mathfunc.h"
#include "md.h"


TINKER_NAMESPACE_BEGIN
void chkpole_acc()
{
   #pragma acc parallel loop independent async deviceptr(x,y,z,zaxis,pole)
   for (int i = 0; i < n; ++i) {
      int polaxe = zaxis[i].polaxe;
      bool check =
         ((polaxe != pole_z_then_x) || (zaxis[i].yaxis) == 0) ? false : true;
      if (check) {
         int k = zaxis[i].yaxis;
         int ia = i;
         int ib = zaxis[i].zaxis;
         int ic = zaxis[i].xaxis;
         int id = INT_ABS(k) - 1;

         // compute the signed parallelpiped volume at chiral site

         real xad = x[ia] - x[id];
         real yad = y[ia] - y[id];
         real zad = z[ia] - z[id];
         real xbd = x[ib] - x[id];
         real ybd = y[ib] - y[id];
         real zbd = z[ib] - z[id];
         real xcd = x[ic] - x[id];
         real ycd = y[ic] - y[id];
         real zcd = z[ic] - z[id];
         real c1 = ybd * zcd - zbd * ycd;
         real c2 = ycd * zad - zcd * yad;
         real c3 = yad * zbd - zad * ybd;
         real vol = xad * c1 + xbd * c2 + xcd * c3;

         // invert atomic multipole components involving the y-axis

         if ((k < 0 && vol > 0) || (k > 0 && vol < 0)) {
            zaxis[i].yaxis = -k;
            // y -> -y
            pole[i][mpl_pme_y] = -pole[i][mpl_pme_y];
            // xy -> -xy
            // yz -> -yz
            pole[i][mpl_pme_xy] = -pole[i][mpl_pme_xy];
            pole[i][mpl_pme_yz] = -pole[i][mpl_pme_yz];
         }
      }
   }
}


#pragma acc routine seq
static void rotsite(int isite, const real (*restrict a)[3],
                    real (*restrict rpole)[10], const real (*restrict pole)[10])
{
   static_assert(mpl_total == 10, "");

   // charge
   rpole[isite][0] = pole[isite][0];
   // dipole
   rpole[isite][1] = 0;
   rpole[isite][2] = 0;
   rpole[isite][3] = 0;
   #pragma acc loop seq collapse(2)
   for (int i = 1; i < 4; ++i) {
      for (int j = 1; j < 4; ++j)
         rpole[isite][i] += pole[isite][j] * a[j - 1][i - 1];
   }
   // quadrupole
   real rp[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
   real mp[3][3];
   mp[0][0] = pole[isite][mpl_pme_xx];
   mp[0][1] = pole[isite][mpl_pme_xy];
   mp[0][2] = pole[isite][mpl_pme_xz];
   mp[1][0] = pole[isite][mpl_pme_yx];
   mp[1][1] = pole[isite][mpl_pme_yy];
   mp[1][2] = pole[isite][mpl_pme_yz];
   mp[2][0] = pole[isite][mpl_pme_zx];
   mp[2][1] = pole[isite][mpl_pme_zy];
   mp[2][2] = pole[isite][mpl_pme_zz];
   #pragma acc loop seq
   for (int i = 0; i < 3; ++i) {
      #pragma acc loop seq
      for (int j = 0; j <= i; ++j) {
         // if (j < i) {
         //  rp[j][i] = rp[i][j];
         // } else {
         #pragma acc loop seq collapse(2)
         for (int k = 0; k < 3; ++k)
            for (int m = 0; m < 3; ++m)
               rp[j][i] += a[k][i] * a[m][j] * mp[m][k];
         // }
      }
   }
   rpole[isite][mpl_pme_xx] = rp[0][0];
   rpole[isite][mpl_pme_xy] = rp[0][1];
   rpole[isite][mpl_pme_xz] = rp[0][2];
   rpole[isite][mpl_pme_yy] = rp[1][1];
   rpole[isite][mpl_pme_yz] = rp[1][2];
   rpole[isite][mpl_pme_zz] = rp[2][2];
}


#pragma acc routine seq
static void rotpole_norm(real* a)
{
   real a1 = REAL_RSQRT(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
   a[0] *= a1;
   a[1] *= a1;
   a[2] *= a1;
}


#pragma acc routine seq
static void rotpole_addto1(real* restrict a, const real* restrict b)
{
   a[0] += b[0];
   a[1] += b[1];
   a[2] += b[2];
}


#pragma acc routine seq
static void rotpole_addto2(real* restrict a, const real* restrict b,
                           const real* restrict c)
{
   a[0] += (b[0] + c[0]);
   a[1] += (b[1] + c[1]);
   a[2] += (b[2] + c[2]);
}


void rotpole_acc()
{
   #pragma acc parallel loop independent async deviceptr(x,y,z,zaxis,rpole,pole)
   for (int i = 0; i < n; ++i) {
      // rotmat routine
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      int iz = zaxis[i].zaxis;
      int ix = zaxis[i].xaxis;
      int iy = INT_ABS(zaxis[i].yaxis) - 1;
      int polaxe = zaxis[i].polaxe;
      // the default identity matrix
      real a[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};


      if (polaxe != pole_none) {
         real* restrict xx = &a[0][0];
         real* restrict yy = &a[1][0];
         real* restrict zz = &a[2][0];


         // STEP 1: PICK Z AND NORM Z
         // pick z
         zz[0] = x[iz] - xi;
         zz[1] = y[iz] - yi;
         zz[2] = z[iz] - zi;
         // norm z
         rotpole_norm(zz);


         // STEP 2: PICK X AND NORM X
         // even if it is not needef for z then x)
         if (polaxe == pole_z_only) {
            // pick x
            int okay = !(REAL_ABS(zz[0]) > 0.866f);
            xx[0] = (okay ? 1 : 0);
            xx[1] = (okay ? 0 : 1);
            xx[2] = 0;
         } else {
            // pick x
            xx[0] = x[ix] - xi;
            xx[1] = y[ix] - yi;
            xx[2] = z[ix] - zi;
            rotpole_norm(xx);
         }


         // STEP 3: PICK Y AND NORM Y
         // only for z biscector and 3 fold
         if (polaxe == pole_z_bisect || polaxe == pole_3_fold) {
            yy[0] = x[iy] - xi;
            yy[1] = y[iy] - yi;
            yy[2] = z[iy] - zi;
            rotpole_norm(yy);
         }


         // STEP 4
         if (polaxe == pole_bisector) {
            rotpole_addto1(zz, xx);
            rotpole_norm(zz);
         } else if (polaxe == pole_z_bisect) {
            rotpole_addto1(xx, yy);
            rotpole_norm(xx);
         } else if (polaxe == pole_3_fold) {
            rotpole_addto2(zz, xx, yy);
            rotpole_norm(zz);
         }


         // STEP 5
         // x -= (x.z) z
         real dotxz = xx[0] * zz[0] + xx[1] * zz[1] + xx[2] * zz[2];
         xx[0] -= dotxz * zz[0];
         xx[1] -= dotxz * zz[1];
         xx[2] -= dotxz * zz[2];
         // norm x
         rotpole_norm(xx);
         // y = z cross x
         a[1][0] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
         a[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
         a[1][2] = a[0][1] * a[2][0] - a[0][0] * a[2][1];
      } // end if (.not. pole_none)


      // rotsite routine
      rotsite(i, a, rpole, pole);
   }
}
TINKER_NAMESPACE_END
