#include "mdpq.h"
#include "tool/io_print.h"
#include <cmath>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/charge.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/dipole.hh>
#include <tinker/detail/moment.hh>
#include <tinker/detail/mpole.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/units.hh>
#include <vector>


namespace tinker {
extern "C" void TINKER_RT(jacobi)(int*, double*, double*, double*);
extern "C" void TINKER_RT(gyrate)(double*);
extern "C" void TINKER_RT(inertia)(int*);
void moments()
{
   moment::netchg = 0;
   moment::netdpl = 0;
   moment::netqpl[0] = 0;
   moment::netqpl[1] = 0;
   moment::netqpl[2] = 0;
   moment::xdpl = 0;
   moment::ydpl = 0;
   moment::zdpl = 0;
   moment::xxqpl = 0;
   moment::xyqpl = 0;
   moment::xzqpl = 0;
   moment::yxqpl = 0;
   moment::yyqpl = 0;
   moment::yzqpl = 0;
   moment::zxqpl = 0;
   moment::zyqpl = 0;
   moment::zzqpl = 0;


   // center of mass
   double weigh = 0, xmid = 0, ymid = 0, zmid = 0;
   std::vector<double> xcm(n), ycm(n), zcm(n);
   for (int i = 0; i < n; ++i) {
      double m = atomid::mass[i];
      weigh += m;
      xmid += atoms::x[i] * m;
      ymid += atoms::y[i] * m;
      zmid += atoms::z[i] * m;
   }
   if (weigh != 0) {
      xmid /= weigh;
      ymid /= weigh;
      zmid /= weigh;
   }
   for (int i = 0; i < n; ++i) {
      xcm[i] = atoms::x[i] - xmid;
      ycm[i] = atoms::y[i] - ymid;
      zcm[i] = atoms::z[i] - zmid;
   }


   // partial charges
   for (int i = 0; i < n and charge::nion > 0; ++i) {
      double c = charge::pchg[i];
      moment::netchg += c;
      moment::xdpl += xcm[i] * c;
      moment::ydpl += ycm[i] * c;
      moment::zdpl += zcm[i] * c;
      moment::xxqpl += xcm[i] * xcm[i] * c;
      moment::xyqpl += xcm[i] * ycm[i] * c;
      moment::xzqpl += xcm[i] * zcm[i] * c;
      moment::yxqpl += ycm[i] * xcm[i] * c;
      moment::yyqpl += ycm[i] * ycm[i] * c;
      moment::yzqpl += ycm[i] * zcm[i] * c;
      moment::zxqpl += zcm[i] * xcm[i] * c;
      moment::zyqpl += zcm[i] * ycm[i] * c;
      moment::zzqpl += zcm[i] * zcm[i] * c;
   }


   // bond dipoles
   for (int i = 0; i < dipole::ndipole; ++i) {
      int j = dipole::idpl[2 * i + 0] - 1;
      int k = dipole::idpl[2 * i + 1] - 1;
      double xi = atoms::x[j] - atoms::x[k];
      double yi = atoms::y[j] - atoms::y[k];
      double zi = atoms::z[j] - atoms::z[k];
      double ri = std::sqrt(xi * xi + yi * yi + zi * zi);
      double xbnd = dipole::bdpl[i] * (xi / ri) / units::debye;
      double ybnd = dipole::bdpl[i] * (yi / ri) / units::debye;
      double zbnd = dipole::bdpl[i] * (zi / ri) / units::debye;
      double xc = atoms::x[j] - xi * dipole::sdpl[i];
      double yc = atoms::y[j] - yi * dipole::sdpl[i];
      double zc = atoms::z[j] - zi * dipole::sdpl[i];
      moment::xdpl += xbnd;
      moment::ydpl += ybnd;
      moment::zdpl += zbnd;
      moment::xxqpl += 2 * xc * xbnd;
      moment::xyqpl += xc * ybnd + yc * xbnd;
      moment::xzqpl += xc * zbnd + zc * xbnd;
      moment::yxqpl += yc * xbnd + xc * ybnd;
      moment::yyqpl += 2 * yc * ybnd;
      moment::yzqpl += yc * zbnd + zc * ybnd;
      moment::zxqpl += zc * xbnd + xc * zbnd;
      moment::zyqpl += zc * ybnd + yc * zbnd;
      moment::zzqpl += 2 * zc * zbnd;
   }


   // atomic multipoles
   for (int i = 0; i < n and mpole::npole > 0; ++i) {
      mpole::rpole[13 * i + 1] += polar::uind[3 * i + 0];
      mpole::rpole[13 * i + 2] += polar::uind[3 * i + 1];
      mpole::rpole[13 * i + 3] += polar::uind[3 * i + 2];


#define RPOLE(j, i) mpole::rpole[13 * (i) + (j)-1]
      moment::netchg += RPOLE(1, i);
      moment::xdpl += xcm[i] * RPOLE(1, i) + RPOLE(2, i);
      moment::ydpl += ycm[i] * RPOLE(1, i) + RPOLE(3, i);
      moment::zdpl += zcm[i] * RPOLE(1, i) + RPOLE(4, i);
      moment::xxqpl += xcm[i] * xcm[i] * RPOLE(1, i) + 2 * xcm[i] * RPOLE(2, i);
      moment::xyqpl += xcm[i] * ycm[i] * RPOLE(1, i) + xcm[i] * RPOLE(3, i) +
         ycm[i] * RPOLE(2, i);
      moment::xzqpl += xcm[i] * zcm[i] * RPOLE(1, i) + xcm[i] * RPOLE(4, i) +
         zcm[i] * RPOLE(2, i);
      moment::yxqpl += +ycm[i] * xcm[i] * RPOLE(1, i) + ycm[i] * RPOLE(2, i) +
         xcm[i] * RPOLE(3, i);
      moment::yyqpl += ycm[i] * ycm[i] * RPOLE(1, i) + 2 * ycm[i] * RPOLE(3, i);
      moment::yzqpl += ycm[i] * zcm[i] * RPOLE(1, i) + ycm[i] * RPOLE(4, i) +
         zcm[i] * RPOLE(3, i);
      moment::zxqpl += zcm[i] * xcm[i] * RPOLE(1, i) + zcm[i] * RPOLE(2, i) +
         xcm[i] * RPOLE(4, i);
      moment::zyqpl += zcm[i] * ycm[i] * RPOLE(1, i) + zcm[i] * RPOLE(3, i) +
         ycm[i] * RPOLE(4, i);
      moment::zzqpl += zcm[i] * zcm[i] * RPOLE(1, i) + 2 * zcm[i] * RPOLE(4, i);
#undef RPOLE
   }


   // convert to traceless quadrupole
   double qave = (moment::xxqpl + moment::yyqpl + moment::zzqpl) / 3;
   moment::xxqpl = 1.5 * (moment::xxqpl - qave);
   moment::xyqpl = 1.5 * moment::xyqpl;
   moment::xzqpl = 1.5 * moment::xzqpl;
   moment::yxqpl = 1.5 * moment::yxqpl;
   moment::yyqpl = 1.5 * (moment::yyqpl - qave);
   moment::yzqpl = 1.5 * moment::yzqpl;
   moment::zxqpl = 1.5 * moment::zxqpl;
   moment::zyqpl = 1.5 * moment::zyqpl;
   moment::zzqpl = 1.5 * (moment::zzqpl - qave);


   // add the atomic quadrupoles
   for (int i = 0; i < n; ++i) {
      moment::xxqpl += 3.0 * mpole::rpole[i * 13 + 4];
      moment::xyqpl += 3.0 * mpole::rpole[i * 13 + 5];
      moment::xzqpl += 3.0 * mpole::rpole[i * 13 + 6];
      moment::yxqpl += 3.0 * mpole::rpole[i * 13 + 7];
      moment::yyqpl += 3.0 * mpole::rpole[i * 13 + 8];
      moment::yzqpl += 3.0 * mpole::rpole[i * 13 + 9];
      moment::zxqpl += 3.0 * mpole::rpole[i * 13 + 10];
      moment::zyqpl += 3.0 * mpole::rpole[i * 13 + 11];
      moment::zzqpl += 3.0 * mpole::rpole[i * 13 + 12];
   }


   // convert dipole to Debye and quadrupole to Buckingham
   moment::xdpl *= units::debye;
   moment::ydpl *= units::debye;
   moment::zdpl *= units::debye;
   moment::xxqpl *= units::debye;
   moment::xyqpl *= units::debye;
   moment::xzqpl *= units::debye;
   moment::yxqpl *= units::debye;
   moment::yyqpl *= units::debye;
   moment::yzqpl *= units::debye;
   moment::zxqpl *= units::debye;
   moment::zyqpl *= units::debye;
   moment::zzqpl *= units::debye;


   // get dipole magnitude and diagonalize quadrupole tensor
   moment::netdpl =
      std::sqrt(moment::xdpl * moment::xdpl + moment::ydpl * moment::ydpl +
                moment::zdpl * moment::zdpl);
   double a[3][3], b[3][3];
   a[0][0] = moment::xxqpl;
   a[1][0] = moment::xyqpl;
   a[2][0] = moment::xzqpl;
   a[0][1] = moment::yxqpl;
   a[1][1] = moment::yyqpl;
   a[2][1] = moment::yzqpl;
   a[0][2] = moment::zxqpl;
   a[1][2] = moment::zyqpl;
   a[2][2] = moment::zzqpl;
   int three = 3;
   TINKER_RT(jacobi)(&three, &a[0][0], moment::netqpl, &b[0][0]);


   auto out = stdout;
   print(out,
         "\n"
         " Total Electric Charge :%12s%13.5lf Electrons\n",
         "", moment::netchg);
   print(out,
         "\n"
         " Dipole Moment Magnitude :%10s%13.3lf Debye\n"
         "\n"
         " Dipole X,Y,Z-Components :%10s%13.3lf%13.3lf%13.3lf\n",
         "", moment::netdpl, "", moment::xdpl, moment::ydpl, moment::zdpl);
   print(out,
         "\n"
         " Quadrupole Moment Tensor :%9s%13.3lf%13.3lf%13.3lf\n"
         "      (Buckinghams)%17s%13.3lf%13.3lf%13.3lf\n"
         "%36s%13.3lf%13.3lf%13.3lf\n",
         "", moment::xxqpl, moment::xyqpl, moment::xzqpl, "", moment::yxqpl,
         moment::yyqpl, moment::yzqpl, "", moment::zxqpl, moment::zyqpl,
         moment::zzqpl

   );
   print(out,
         "\n"
         " Principal Axes Quadrupole :%8s%13.3lf%13.3lf%13.3lf\n",
         "", moment::netqpl[0], moment::netqpl[1], moment::netqpl[2]);


   if (chgpot::dielec != 1) {
      print(out,
            "\n"
            " Dielectric Constant :%14s%13.3lf\n",
            "", chgpot::dielec);
      print(out, " Effective Total Charge :%11s%13.5lf Electrons\n", "",
            moment::netchg / std::sqrt(chgpot::dielec));
      print(out, " Effective Dipole Moment :%10s%13.3lf Debye\n", "",
            moment::netdpl / std::sqrt(chgpot::dielec));
   }


   // radius of gyration and moments of inertia
   double rg;
   TINKER_RT(gyrate)(&rg);
   print(out,
         "\n"
         " Radius of Gyration :%15s%13.3lf Angstroms\n",
         "", rg);
   int one = 1;
   TINKER_RT(inertia)(&one);
}
}
