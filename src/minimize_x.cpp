#include "energy.h"
#include "fcxx.h"
#include "io_fort_str.h"
#include "io_print.h"
#include "nblist.h"
#include "tinker_rt.h"
#include <tinker/detail/files.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/scales.hh>


TINKER_NAMESPACE_BEGIN
namespace {
std::vector<double> grx, gry, grz;
}


double minimiz1(double* xx, double* g);


void minimize_set_xx(int n, double* xx, const double* scale);
void minimize_set_xyz(int n, const double* xx, const double* scale);
void minimize_set_xx_by_pos_acc(int, double*, const double*);
void minimize_set_pos_acc(int, const double*, const double*);


void x_minimize(int argc, char** argv)
{
   TINKER_RT(initial)();
   TINKER_RT(getxyz)();
   TINKER_RT(mechanic)();


   // perform the setup functions needed for optimization
   t_optinit();


   // get termination criterion as RMS gradient per atom
   int exist = false;
   char string[240];
   double grdmin = -1.0;
   nextarg(string, exist);
   if (exist) {
      read_string(grdmin, string);
   }
   std::string prompt = "\n"
                        " Enter RMS Gradient per Atom Criterion [0.01] :  ";
   read_stream(grdmin, prompt, 0.01, [](double val) { return val < 0; });


   // write out a copy of coordinates for later update
   int imin = t_freeunit();
   const int leng = files::leng;
   std::string minfile = fstr_view(files::filename)(1, leng).trim() + ".xyz";
   minfile = t_version(minfile, "new");
   t_open(imin, minfile.c_str(), "new");
   t_prtxyz(imin);
   t_close(imin);
   fstr_view outview = files::outfile;
   outview = minfile;


   int flags = calc::xyz + calc::mass;
   flags += (calc::energy + calc::grad);
   rc_flag = flags;
   initialize();


   // perform dynamic allocation of some global arrays
   if (!t_allocated(scales::scale))
      t_allocate_d1(&scales::scale, 3 * n);


   // set scaling parameter for function and derivative values;
   // use square root of median eigenvalue of typical Hessian
   scales::set_scale = true;
   int nvar = 0;
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < 3; ++j) {
         scales::scale[nvar] = 12;
         ++nvar;
      }
   }


   // perform dynamic allocation of some local arrays
   std::vector<double> xxvec(3 * n);
   double* xx = xxvec.data();
   grx.resize(n);
   gry.resize(n);
   grz.resize(n);


   // convert atomic coordinates to optimization parameters
   minimize_set_xx(n, xx, scales::scale);


   // make the call to the optimization routine
   t_lbfgs(3 * n, xx, grdmin, (void*)minimiz1);


   // convert optimization parameters to atomic coordinates
   minimize_set_xyz(n, xx, scales::scale);


   // compute the final function and RMS gradient values
   double minimum;
   energy(calc::energy + calc::grad);
   copy_energy(calc::energy, &minimum);
   copy_gradient(calc::grad, grx.data(), gry.data(), grz.data());
   double gnorm = 0;
   for (int i = 0; i < n; ++i) {
      gnorm += grx[i] * grx[i] + gry[i] * gry[i] + grz[i] * grz[i];
   }
   gnorm = std::sqrt(gnorm);
   double grms = gnorm / std::sqrt(n);


   // perform deallocation of some local arrays
   xxvec.clear();
   grx.clear();
   gry.clear();
   grz.clear();


   // write out the final function and gradient values
   auto o = stdout;
   const int d1 = inform::digits;
   int l1 = 12 + d1;
   auto fstr = "{0:{1}.{2}f}"_s;
   auto estr = fstr;
   if (grms <= std::pow(10, -d1))
      estr = "{0:{1}.{2}E}"_s;
   print(o, "\n Final Function Value :  "_s + fstr, minimum, l1, d1);
   print(o, "\n Final RMS Gradient :    "_s + estr, grms, l1, d1);
   print(o, "\n Final Gradient Norm :   "_s + estr, gnorm, l1, d1);
   print(o, "\n");


   // write the final coordinates into a file
   bounds();
   imin = t_freeunit();
   t_open(imin, minfile, "old");
   t_rewind(imin);
   t_prtxyz(imin);
   t_close(imin);


   finish();
   TINKER_RT(final)();
}


double minimiz1(double* xx, double* g)
{
   // convert optimization parameters to atomic coordinates
   minimize_set_xyz(n, xx, scales::scale);


   // compute and store the energy and gradient
   energy(calc::energy + calc::grad);
   energy_prec eout;
   copy_energy(calc::energy, &eout);
   copy_gradient(calc::grad, grx.data(), gry.data(), grz.data());


   // convert coordinates and gradient to optimization parameters
   minimize_set_xx(n, xx, scales::scale);
   for (int i = 0; i < n; ++i) {
      int ii = 3 * i;
      g[ii + 0] = grx[i] / scales::scale[ii + 0];
      g[ii + 1] = gry[i] / scales::scale[ii + 1];
      g[ii + 2] = grz[i] / scales::scale[ii + 2];
   }


   return eout;
}


void minimize_set_xx(int n, double* xx, const double* scale)
{
   minimize_set_xx_by_pos_acc(n, xx, scale);
}


void minimize_set_xyz(int n, const double* xx, const double* scale)
{
   minimize_set_pos_acc(n, xx, scale);
   copy_pos_to_xyz();
   refresh_neighbors();
}
TINKER_NAMESPACE_END
