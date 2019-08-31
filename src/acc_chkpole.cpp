#include "acc_mathfunc.h"
#include "e_mpole.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
void chkpole() {
  auto* pole = pole_vec.data();
  auto* rpole = rpole_vec.data();
  #pragma acc data deviceptr(x,y,z,zaxis,pole)
  #pragma acc parallel loop
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
TINKER_NAMESPACE_END
