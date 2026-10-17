#include "ff/box.h"
#include "ff/energy.h"
#include "ff/nblist.h"
#include "ff/rwcrd.h"
#include "math/random.h"
#include "tool/argkey.h"
#include "tool/error.h"
#include "tool/iofortstr.h"
#include "tool/ioprint.h"
#include "tool/ioread.h"
#include <tinker/detail/files.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/keys.hh>
#include <tinker/detail/output.hh>
#include <tinker/detail/titles.hh>
#include <tinker/detail/units.hh>
#include <tinker/routines.h>

#include "tinker9.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <tuple>

namespace tinker {
static constexpr int MAX_NCHAR = 240;

static void xBarMake()
{
   rc_flag = calc::xyz | calc::mass | calc::energy;

   auto out = stdout;
   auto invalid_temperature = [](double t) { return t <= 0; };
   auto invalid_recompute_answer = [](char c) { return c != 'y' and c != 'Y' and c != 'n' and c != 'N'; };

   int lenga, lengb;
   int ltitlea, ltitleb;
   int nkey0, nkey1;
   double tempa, tempb;
   char filea[MAX_NCHAR], fileb[MAX_NCHAR];
   char titlea[MAX_NCHAR], titleb[MAX_NCHAR];
   std::vector<std::array<char, MAX_NCHAR>> keys0, keys1;
   char string[MAX_NCHAR];
   std::string str;

   // get trajectory A and setup mechanics calculation
   int iarc, first;
   tinker_f_getarc(&iarc);
   tinker_f_close(&iarc);
   tinker_f_mechanic();
   mechanic2();

   // store the filename for trajectory A
   std::memcpy(filea, files::filename, MAX_NCHAR);
   lenga = files::leng;
   std::memcpy(titlea, titles::title, MAX_NCHAR);
   ltitlea = titles::ltitle;

   // store the keyword values for state 0
   nkey0 = keys::nkey;
   keys0.resize(nkey0);
   for (int i = 0; i < nkey0; ++i)
      std::memcpy(keys0[i].data(), keys::keyline[i], MAX_NCHAR);

   // find the original temperature value for trajectory A
   bool exist = false;
   tempa = -1;
   nextarg(string, exist);
   if (exist)
      ioReadString(tempa, string);
   ioReadStream(tempa,
      "\n"
      " Enter Trajectory A Temperature in Degrees K [298] :  ",
      298.0, invalid_temperature);

   // get trajectory B and setup mechanics calculation
   tinker_f_getarc(&iarc);
   tinker_f_close(&iarc);
   tinker_f_mechanic();
   mechanic2();
   inform::silent = 1;

   // store the filename for trajectory B
   std::memcpy(fileb, files::filename, MAX_NCHAR);
   lengb = files::leng;
   std::memcpy(titleb, titles::title, MAX_NCHAR);
   ltitleb = titles::ltitle;

   // store the keyword values for state 1
   nkey1 = keys::nkey;
   keys1.resize(nkey1);
   for (int i = 0; i < nkey1; ++i)
      std::memcpy(keys1[i].data(), keys::keyline[i], MAX_NCHAR);

   // find the original temperature value for trajectory B
   tempb = -1;
   nextarg(string, exist);
   if (exist)
      ioReadString(tempb, string);
   ioReadStream(tempb,
      "\n"
      " Enter Trajectory B Temperature in Degrees K [298] :  ",
      298.0, invalid_temperature);

   // decide whether to use energies from trajectory log files
   bool recompute = true;
   char answer = ' ';
   nextarg(string, exist);
   if (exist)
      ioReadString(answer, string);
   ioReadStream(answer,
      "\n"
      " Obtain Energies from Trajectory Logs if Available [N] :  ",
      'N', invalid_recompute_answer);
   if (answer == 'y' or answer == 'Y')
      recompute = false;

   //====================================================================//

   int done;
   const char* log_fmt = " Current Potential %lf\n";
   const char* process_fmt = "       Completed%8d Coordinate Frames\n";
   // the buffer must hold only the blank-padded base name, or suffix finds the
   // extension of whatever it held before
   std::memset(string, ' ', MAX_NCHAR);
   std::memcpy(string, filea, lenga);
   tinker_f_suffix({string, MAX_NCHAR}, {const_cast<char*>("bar"), 3}, {const_cast<char*>("new"), 3});
   std::string barfile = FstrView(string).trim();
   std::vector<double> ua0, ua1, ub0, ub1, vola, volb;

   // check for log with energies of trajectory A in state 0
   if (not recompute) {
      std::memset(string, ' ', MAX_NCHAR);
      std::memcpy(string, filea, lenga);
      tinker_f_suffix({string, MAX_NCHAR}, {const_cast<char*>("log"), 3}, {const_cast<char*>("old"), 3});
      str = FstrView(string).trim();
      std::ifstream a_log(str);
      while (std::getline(a_log, str)) {
         double val;
         int count = std::sscanf(str.c_str(), log_fmt, &val);
         if (count == 1)
            ua0.push_back(val);
      }
      a_log.close();
   }

   // reopen the file corresponding to trajectory A
   print(out,
      "\n"
      " Initial Processing for Trajectory A :\n"
      "\n");

   std::memcpy(string, filea, MAX_NCHAR);
   tinker_f_suffix({string, MAX_NCHAR}, {const_cast<char*>("arc"), 3}, {const_cast<char*>("old"), 3});
   iarc = tinker_f_freeunit();
   str = FstrView(string).trim();
   if (output::archive)
      tinker_f_open(&iarc, str, "old");
   else if (output::binary)
      tinker_f_open(&iarc, str, "unformatted", "old");

   if (ua0.size() == 0) {
      // reset trajectory A using the parameters for state 0
      tinker_f_rewind(&iarc);
      first = 1;
      tinker_f_readcart(&iarc, &first);
      keys::nkey = nkey0;
      for (int i = 0; i < keys::nkey; ++i)
         std::memcpy(keys::keyline[i], keys0[i].data(), MAX_NCHAR);
      tinker_f_mechanic();
      mechanic2();

      // find potential energies for trajectory A in state 0
      initialize();
      done = 0;
      auto a_arc = CrdReader(str);
      do {
         done = a_arc.readCurrent();
         nblistRefresh();
         energy(calc::energy);
         ua0.push_back(esum);
         int i = ua0.size();
         if (i % 100 == 0) {
            print(out, process_fmt, i);
            std::fflush(out);
         }
      } while (not done);
      finish();
   } else {
      int ii = ua0.size();
      for (int i = 1; i <= ii; ++i) {
         if (i % 100 == 0)
            print(out, process_fmt, i);
      }
      std::fflush(out);
   }

   // reset trajectory A using the parameters for state 1
   tinker_f_rewind(&iarc);
   first = 1;
   tinker_f_readcart(&iarc, &first);
   keys::nkey = nkey1;
   for (int i = 0; i < keys::nkey; ++i)
      std::memcpy(keys::keyline[i], keys1[i].data(), MAX_NCHAR);
   tinker_f_mechanic();
   mechanic2();

   // find potential energies for trajectory A in state 1
   if (inform::verbose)
      print(out,
         "\n"
         " Potential Energy Values for Trajectory A :"
         "\n\n"
         "       Frame         State 0         State 1            Delta\n"
         "\n");
   initialize();
   done = 0;
   auto a_ar2 = CrdReader(str);
   do {
      done = a_ar2.readCurrent();
      nblistRefresh();
      energy(calc::energy);
      // a system without a periodic box has no volume
      double vol = box_shape == BoxShape::UNBOUND ? 0 : boxVolume();
      ua1.push_back(esum);
      vola.push_back(vol);
      if (inform::verbose and ua1.size() <= ua0.size()) {
         int i = ua1.size() - 1;
         print(out, "%11d  %16.4lf%16.4lf%16.4lf\n", i + 1, ua0[i], ua1[i], ua1[i] - ua0[i]);
      }
   } while (not done);
   finish();

   // save potential energies and volumes for trajectory A
   tinker_f_close(&iarc);
   FILE* ibar = std::fopen(barfile.c_str(), "w");
   int nfrma = std::min(ua0.size(), ua1.size());
   str = FstrView(titlea)(1, ltitlea).trim();
   print(ibar, "%8d%10.2lf  %s\n", nfrma, tempa, str);
   for (int i = 0; i < nfrma; ++i) {
      if (vola[i] == 0)
         print(ibar, "%8d  %18.4lf%18.4lf\n", i + 1, ua0[i], ua1[i]);
      else
         print(ibar, "%8d  %18.4lf%18.4lf%18.4lf\n", i + 1, ua0[i], ua1[i], vola[i]);
   }
   std::fflush(ibar);

   //====================================================================//

   // check for log with energies of trajectory B in state 1
   if (not recompute) {
      std::memset(string, ' ', MAX_NCHAR);
      std::memcpy(string, fileb, lengb);
      tinker_f_suffix({string, MAX_NCHAR}, {const_cast<char*>("log"), 3}, {const_cast<char*>("old"), 3});
      str = FstrView(string).trim();
      std::ifstream b_log(str);
      while (std::getline(b_log, str)) {
         double val;
         int count = std::sscanf(str.c_str(), log_fmt, &val);
         if (count == 1)
            ub1.push_back(val);
      }
      b_log.close();
   }

   // reopen the file corresponding to trajectory B
   print(out,
      "\n"
      " Initial Processing for Trajectory B :\n"
      "\n");

   std::memcpy(string, fileb, MAX_NCHAR);
   tinker_f_suffix({string, MAX_NCHAR}, {const_cast<char*>("arc"), 3}, {const_cast<char*>("old"), 3});
   iarc = tinker_f_freeunit();
   str = FstrView(string).trim();
   if (output::archive)
      tinker_f_open(&iarc, str, "old");
   else if (output::binary)
      tinker_f_open(&iarc, str, "unformatted", "old");

   if (ub1.size() == 0) {
      // reset trajectory B using the parameters for state 1
      tinker_f_rewind(&iarc);
      first = 1;
      tinker_f_readcart(&iarc, &first);
      keys::nkey = nkey1;
      for (int i = 0; i < keys::nkey; ++i)
         std::memcpy(keys::keyline[i], keys1[i].data(), MAX_NCHAR);
      tinker_f_mechanic();
      mechanic2();

      // find potential energies for trajectory B in state 1
      initialize();
      done = 0;
      auto b_arc = CrdReader(str);
      do {
         done = b_arc.readCurrent();
         nblistRefresh();
         energy(calc::energy);
         ub1.push_back(esum);
         int i = ub1.size();
         if (i % 100 == 0) {
            print(out, process_fmt, i);
            std::fflush(out);
         }
      } while (not done);
      finish();
   } else {
      int ii = ub1.size();
      for (int i = 1; i <= ii; ++i) {
         if (i % 100 == 0)
            print(out, process_fmt, i);
      }
      std::fflush(out);
   }

   // reset trajectory B using the parameters for state 0
   tinker_f_rewind(&iarc);
   first = 1;
   tinker_f_readcart(&iarc, &first);
   keys::nkey = nkey0;
   for (int i = 0; i < keys::nkey; ++i)
      std::memcpy(keys::keyline[i], keys0[i].data(), MAX_NCHAR);
   tinker_f_mechanic();
   mechanic2();

   // find potential energies for trajectory B in state 0
   if (inform::verbose)
      print(out,
         "\n"
         " Potential Energy Values for Trajectory B :"
         "\n\n"
         "       Frame         State 0         State 1            Delta\n"
         "\n");
   initialize();
   done = 0;
   auto b_ar2 = CrdReader(str);
   do {
      done = b_ar2.readCurrent();
      nblistRefresh();
      energy(calc::energy);
      double vol = box_shape == BoxShape::UNBOUND ? 0 : boxVolume();
      ub0.push_back(esum);
      volb.push_back(vol);
      if (inform::verbose and ub0.size() <= ub1.size()) {
         int i = ub0.size() - 1;
         print(out, "%11d  %16.4lf%16.4lf%16.4lf\n", i + 1, ub0[i], ub1[i], ub0[i] - ub1[i]);
      }
   } while (not done);
   finish();

   // save potential energies and volumes for trajectory B
   tinker_f_close(&iarc);
   int nfrmb = std::min(ub0.size(), ub1.size());
   str = FstrView(titleb)(1, ltitleb).trim();
   print(ibar, "%8d%10.2lf  %s\n", nfrmb, tempb, str);
   for (int i = 0; i < nfrmb; ++i) {
      if (volb[i] == 0)
         print(ibar, "%8d  %18.4lf%18.4lf\n", i + 1, ub0[i], ub1[i]);
      else
         print(ibar, "%8d  %18.4lf%18.4lf%18.4lf\n", i + 1, ub0[i], ub1[i], volb[i]);
   }
   std::fclose(ibar);
   print(out,
      "\n"
      " Potential Energy Values Written To :  %s\n",
      barfile);
}

namespace {
/// Potential energies and volumes of the frames selected from both trajectories.
struct BarData
{
   std::vector<double> ua0, ua1; ///< Energies of trajectory A frames in states 0 and 1.
   std::vector<double> ub0, ub1; ///< Energies of trajectory B frames in states 0 and 1.
   std::vector<double> vola;     ///< Volumes of trajectory A frames, or 0 if not used.
   std::vector<double> volb;     ///< Volumes of trajectory B frames, or 0 if not used.
};

/// Fills idx with n frame indices drawn with replacement.
void barDraw(int n, std::vector<int>& idx)
{
   idx.resize(n);
   for (int i = 0; i < n; ++i)
      idx[i] = std::min((int)(n * random<double>()), n - 1);
}

/// Finds the free energy and enthalpy change from the sampled state to the
/// other state via thermodynamic perturbation over a list of frames, which may
/// repeat, along with the average energy and volume of those frames. The
/// exponentials are shifted by their largest exponent so the sum cannot
/// overflow, and the enthalpy includes the PV term.
void barFep(const std::vector<int>& idx, const std::vector<double>& uo, const std::vector<double>& ut,
   const std::vector<double>& vol, double rt, double pres, double& dg, double& dh, double& uave, double& vave)
{
   const int n = idx.size();
   double xmax = 0;
   for (int i = 0; i < n; ++i) {
      int j = idx[i];
      double x = (uo[j] - ut[j]) / rt;
      if (i == 0 or x > xmax)
         xmax = x;
   }

   double wsum = 0, hwsum = 0, usum = 0, vsum = 0;
   for (int i = 0; i < n; ++i) {
      int j = idx[i];
      double w = std::exp((uo[j] - ut[j]) / rt - xmax);
      wsum += w;
      hwsum += (ut[j] + pres * vol[j]) * w;
      usum += uo[j];
      vsum += vol[j];
   }
   uave = usum / n;
   vave = vsum / n;
   dg = -rt * (xmax + std::log(wsum / n));
   dh = hwsum / wsum - (uave + pres * vave);
}

/// Solves the Bennett acceptance ratio equation for the free energy change from
/// state 0 to state 1 by iteration over lists of frames sampled in each state,
/// starting from cinit; also finds the analytical error of Bennett. Records the
/// value at each iteration in hist if it is not null, and returns false if the
/// iteration does not converge.
bool barIter(const std::vector<int>& ia, const std::vector<int>& ib, const BarData& d, double rta, double rtb,
   double cinit, int maxiter, double& cnew, double& stdev, std::vector<double>* hist)
{
   const double eps = 0.0001;
   const double rt = 0.5 * (rta + rtb);
   const int na = ia.size(), nb = ib.size();
   const double frma = na, frmb = nb;
   const double rfrm = frma / frmb;

   cnew = cinit;
   stdev = 0;
   if (hist)
      hist->clear();
   for (int iter = 0;; ++iter) {
      double cold = cnew;
      double top = 0, top2 = 0;
      for (int i = 0; i < nb; ++i) {
         int j = ib[i];
         double fterm = 1.0 / (1.0 + std::exp((d.ub0[j] - d.ub1[j] + cold) / rtb));
         top += fterm;
         top2 += fterm * fterm;
      }
      double bot = 0, bot2 = 0;
      for (int i = 0; i < na; ++i) {
         int j = ia[i];
         double fterm = 1.0 / (1.0 + std::exp((d.ua1[j] - d.ua0[j] - cold) / rta));
         bot += fterm;
         bot2 += fterm * fterm;
      }
      cnew = rt * std::log(rfrm * top / bot) + cold;
      stdev = rt * std::sqrt(std::max(0.0, (bot2 - bot * bot / frma) / (bot * bot) + (top2 - top * top / frmb) / (top * top)));
      if (hist)
         hist->push_back(cnew);
      if (std::fabs(cnew - cold) < eps)
         return true;
      else if (iter >= maxiter)
         return false;
      else if (not(std::fabs(cnew) < 1.0e100))
         return false;
   }
}

/// Finds the enthalpy change from state 0 to state 1 via the BAR estimator of
/// Wyczalkowski et al. over lists of frames sampled in each state, given the BAR
/// free energy change cbar for the same frames. The state energies include the
/// PV term, which cancels in the energy differences.
double barHbar(const std::vector<int>& ia, const std::vector<int>& ib, const BarData& d, double rta, double rtb,
   double pres, double cbar)
{
   // find the terms from the state 0 frames
   double fsum = 0, fvsum = 0, fbvsum = 0, vsum = 0, fbsum0 = 0;
   for (int j : ia) {
      double fore = 1.0 / (1.0 + std::exp((d.ua1[j] - d.ua0[j] - cbar) / rta));
      double back = 1.0 / (1.0 + std::exp((d.ua0[j] - d.ua1[j] + cbar) / rta));
      double h = d.ua0[j] + pres * d.vola[j];
      fsum += fore;
      fvsum += fore * h;
      fbvsum += fore * back * (d.ua1[j] - d.ua0[j]);
      vsum += h;
      fbsum0 += fore * back;
   }
   double alpha0 = fvsum - fsum * (vsum / ia.size()) + fbvsum;

   // find the terms from the state 1 frames
   double bsum = 0, bvsum = 0, fbsum1 = 0;
   fbvsum = 0, vsum = 0;
   for (int j : ib) {
      double fore = 1.0 / (1.0 + std::exp((d.ub1[j] - d.ub0[j] - cbar) / rtb));
      double back = 1.0 / (1.0 + std::exp((d.ub0[j] - d.ub1[j] + cbar) / rtb));
      double h = d.ub1[j] + pres * d.volb[j];
      bsum += back;
      bvsum += back * h;
      fbvsum += fore * back * (d.ub1[j] - d.ub0[j]);
      vsum += h;
      fbsum1 += fore * back;
   }
   double alpha1 = bvsum - bsum * (vsum / ib.size()) - fbvsum;
   return (alpha0 - alpha1) / (fbsum0 + fbsum1);
}

/// Computes and prints the free energy, enthalpy and entropy change from state 0
/// to state 1 via the FEP and BAR methods, with the PV term at patm atmospheres.
/// Each value is computed once from all of the frames, and the BAR enthalpy and
/// entropy use the BAR iteration free energy. The BAR iteration error is the
/// analytical estimate of Bennett; the other errors are standard deviations over
/// nbst bootstrap trials that apply the same estimators to frames resampled with
/// replacement. The BAR bootstrap value is the mean over converged trials only,
/// and BAR errors are reported only if at least two trials converge. Throws if
/// the BAR iteration over all of the frames does not converge.
void barPair(const BarData& d, double tempa, double tempb, double patm, int nbst)
{
   auto out = stdout;
   const int nfrma = d.ua0.size(), nfrmb = d.ub0.size();
   const double rta = units::gasconst * tempa;
   const double rtb = units::gasconst * tempb;
   const double pres = patm / units::prescon;
   const double bst = nbst;
   const double ratio = bst / (bst - 1);
   const int maxiter = 100;

   // find the FEP free energies, enthalpies and averages from all of the
   // frames, where the enthalpies include the PV term
   std::vector<int> ia(nfrma), ib(nfrmb);
   for (int i = 0; i < nfrma; ++i)
      ia[i] = i;
   for (int i = 0; i < nfrmb; ++i)
      ib[i] = i;
   double cfore, hfore, uave0, vavea;
   double cback, hback, uave1, vaveb;
   barFep(ia, d.ua0, d.ua1, d.vola, rta, pres, cfore, hfore, uave0, vavea);
   barFep(ib, d.ub1, d.ub0, d.volb, rtb, pres, cback, hback, uave1, vaveb);
   cback = -cback;
   hback = -hback;
   double epv = pres * (vaveb - vavea);
   double hdir = uave1 - uave0 + epv;

   // find the BAR free energy from all of the frames, starting from the mean
   // of the FEP values, then the BAR enthalpy at that value
   double cinit = 0.5 * (cfore + cback);
   if (not(std::fabs(cinit) < 1.0e100))
      cinit = 0;
   double cbar, stdev;
   std::vector<double> hist;
   bool conv = barIter(ia, ib, d, rta, rtb, cinit, maxiter, cbar, stdev, &hist);
   double hbr = 0;
   if (conv)
      hbr = barHbar(ia, ib, d, rta, rtb, pres, cbar);

   // apply the same estimators to bootstrap resamples of the frames, with each
   // BAR trial started from the converged BAR value and dropped if it fails to
   // converge
   double val[8], sd[8], vsum[8] = {0}, vsum2[8] = {0};
   int nconv = 0;
   double csum = 0, csum2 = 0, hsum = 0, hsum2 = 0;
   for (int k = 0; k < nbst; ++k) {
      barDraw(nfrma, ia);
      barDraw(nfrmb, ib);
      double vtrla, vtrlb;
      barFep(ia, d.ua0, d.ua1, d.vola, rta, pres, val[0], val[6], val[2], vtrla);
      barFep(ib, d.ub1, d.ub0, d.volb, rtb, pres, val[1], val[7], val[3], vtrlb);
      val[1] = -val[1];
      val[7] = -val[7];
      val[4] = pres * (vtrlb - vtrla);
      val[5] = val[3] - val[2] + val[4];
      for (int m = 0; m < 8; ++m) {
         vsum[m] += val[m];
         vsum2[m] += val[m] * val[m];
      }
      if (conv) {
         double ctrial, etrial;
         if (barIter(ia, ib, d, rta, rtb, cbar, maxiter, ctrial, etrial, nullptr)) {
            ++nconv;
            csum += ctrial;
            csum2 += ctrial * ctrial;
            double htrial = barHbar(ia, ib, d, rta, rtb, pres, ctrial);
            hsum += htrial;
            hsum2 += htrial * htrial;
         }
      }
   }
   for (int m = 0; m < 8; ++m) {
      double ave = vsum[m] / bst;
      sd[m] = std::sqrt(std::max(0.0, ratio * (vsum2[m] / bst - ave * ave)));
   }

   // report the free energy change via thermodynamic perturbation
   print(out,
      "\n"
      " Free Energy Difference via FEP Method :\n"
      "\n");
   print(out, " Free Energy via Forward FEP         %12.4lf +/-%9.4lf Kcal/mol\n", cfore, sd[0]);
   print(out, " Free Energy via Backward FEP        %12.4lf +/-%9.4lf Kcal/mol\n", cback, sd[1]);

   // report the free energy change via the BAR iteration
   print(out,
      "\n"
      " Free Energy Difference via BAR Method :\n"
      "\n");
   const int niter = (int)hist.size() - 1;
   for (int i = 0; i <= niter; ++i)
      print(out, " BAR Iteration%4d%19s%12.4lf Kcal/mol\n", i, "", hist[i]);
   if (conv and niter == 0)
      print(out, " BAR Free Energy Estimate        %12.4lf +/-%9.4lf Kcal/mol\n", cbar, stdev);
   else if (conv)
      print(out,
         "\n"
         " Free Energy via BAR Estimate        %12.4lf +/-%9.4lf Kcal/mol\n",
         cbar, stdev);

   // stop if the BAR iteration fails
   if (not conv) {
      print(out,
         "\n"
         " BAR Free Energy Estimate not Converged after%4d Iterations\n",
         niter);
      std::fflush(out);
      TINKER_THROW("BAR  --  Free Energy Estimate not Converged");
   }

   // report the mean and spread of the BAR bootstrap trials
   const double rconv = nconv;
   double stdhbr = 0;
   if (nconv >= 2) {
      double cave = csum / rconv;
      double cstd = std::sqrt(std::max(0.0, rconv / (rconv - 1) * (csum2 / rconv - cave * cave)));
      double have = hsum / rconv;
      stdhbr = std::sqrt(std::max(0.0, rconv / (rconv - 1) * (hsum2 / rconv - have * have)));
      print(out, " Free Energy via BAR Bootstrap       %12.4lf +/-%9.4lf Kcal/mol\n", cave, cstd);
   } else {
      print(out, " Free Energy via BAR Bootstrap             not Available\n");
   }
   if (nconv < nbst)
      print(out, " BAR Bootstrap :%8d of%8d Trials did not Converge\n", nbst - nconv, nbst);

   // report the enthalpy directly via average potential energy
   print(out,
      "\n"
      " Enthalpy from Potential Energy Averages :\n"
      "\n");
   print(out, " Average Energy for State 0          %12.4lf +/-%9.4lf Kcal/mol\n", uave0, sd[2]);
   print(out, " Average Energy for State 1          %12.4lf +/-%9.4lf Kcal/mol\n", uave1, sd[3]);
   if (epv != 0)
      print(out, " PdV Work Term at%9.2lf Atm       %12.4lf +/-%9.4lf Kcal/mol\n", patm, epv, sd[4]);
   print(out, " Enthalpy via Direct Estimate        %12.4lf +/-%9.4lf Kcal/mol\n", hdir, sd[5]);

   // report the enthalpy and entropy via thermodynamic perturbation
   double sfore = (hfore - cfore) / tempa;
   double sback = (hback - cback) / tempb;
   print(out,
      "\n"
      " Enthalpy and Entropy via FEP Method :\n"
      "\n");
   print(out, " Enthalpy via Forward FEP            %12.4lf +/-%9.4lf Kcal/mol\n", hfore, sd[6]);
   print(out, " Entropy via Forward FEP             %12.6lf Kcal/mol/K\n", sfore);
   print(out, " Forward FEP -T*dS Value             %12.4lf Kcal/mol\n", -tempa * sfore);
   print(out,
      "\n"
      " Enthalpy via Backward FEP           %12.4lf +/-%9.4lf Kcal/mol\n",
      hback, sd[7]);
   print(out, " Entropy via Backward FEP            %12.6lf Kcal/mol/K\n", sback);
   print(out, " Backward FEP -T*dS Value            %12.4lf Kcal/mol\n", -tempb * sback);

   // the BAR enthalpy and entropy need converged bootstrap trials
   print(out,
      "\n"
      " Enthalpy and Entropy via BAR Method :\n"
      "\n");
   if (nconv == 0) {
      print(out, " BAR Enthalpy and Entropy not Computed\n");
      return;
   }

   // report the enthalpy and entropy via the BAR method
   double tsbar = hbr - cbar;
   double sbar = tsbar / (0.5 * (tempa + tempb));
   print(out, " Enthalpy via BAR Estimate           %12.4lf +/-%9.4lf Kcal/mol\n", hbr, stdhbr);
   print(out, " Entropy via BAR Estimate            %12.6lf Kcal/mol/K\n", sbar);
   print(out, " BAR Estimate -T*dS Value            %12.4lf Kcal/mol\n", -tsbar);
}
}

// Mode 2 finds thermodynamic values from a BAR file via barPair. Volumes in the
// file enter only the PV terms of the enthalpies, at a pressure taken from the
// command line or asked for, and both trajectories must be at one temperature.
static void xBarCalc()
{
   auto out = stdout;
   bool exist = false;
   char string[MAX_NCHAR];
   std::string str;
   // the BAR file a name refers to, with the "bar" extension if it has none
   auto bar_path = [&](const std::string& name) {
      std::memset(string, ' ', MAX_NCHAR);
      std::memcpy(string, name.data(), std::min(name.length(), (size_t)MAX_NCHAR));
      tinker_f_basefile({string, MAX_NCHAR});
      tinker_f_suffix({string, MAX_NCHAR}, {const_cast<char*>("bar"), 3}, {const_cast<char*>("old"), 3});
      return std::string(FstrView(string).trim());
   };
   auto invalid_barfile = [&](const std::string& name) { return not std::ifstream(bar_path(name)); };

   // ask the user for file with potential energies and volumes
   nextarg(string, exist);
   if (exist)
      str = FstrView(string).trim();
   ioReadStream(str,
      "\n"
      " Enter Potential Energy BAR File Name :  ",
      std::string(""), invalid_barfile);
   std::string barfile = bar_path(str);

   // a frame selection is invalid only if a value is negative; a zero or
   // missing value selects the first frame, the last frame or a step of one
   auto invalid_frames = [](const std::string& s) {
      int v[3] = {0, 0, 0};
      std::istringstream iss(s);
      for (int k = 0; k < 3 and (iss >> v[k]); ++k) {}
      return v[0] < 0 or v[1] < 0 or v[2] < 0;
   };
   auto get_frames = [&](const char* prompt, int& start, int& stop, int& step) {
      start = 0, stop = 0, step = 0;
      bool query = true;
      nextarg(string, exist);
      if (exist and std::sscanf(string, "%d", &start) == 1) {
         query = false;
         nextarg(string, exist);
         if (exist and std::sscanf(string, "%d", &stop) == 1) {
            nextarg(string, exist);
            if (exist)
               std::sscanf(string, "%d", &step);
         }
      }
      if (query)
         str = "-1";
      else
         str = std::to_string(start) + " " + std::to_string(stop) + " " + std::to_string(step);
      ioReadStream(str, prompt, std::string(""), invalid_frames);
      start = 0, stop = 0, step = 0;
      std::istringstream iss(str);
      iss >> start >> stop >> step;
      if (start == 0)
         start = 1;
      if (stop == 0)
         stop = std::numeric_limits<int>::max();
      if (step == 0)
         step = 1;
   };

   // set beginning and ending frame for trajectory A
   int starta, stopa, stepa;
   get_frames("\n"
              " First & Last Frame and Step Increment for Trajectory A :  ",
      starta, stopa, stepa);

   // set beginning and ending frame for trajectory B
   int startb, stopb, stepb;
   get_frames("\n"
              " First & Last Frame and Step Increment for Trajectory B :  ",
      startb, stopb, stepb);

   // Reads one trajectory block of the BAR file, a header with the frame
   // count, temperature and title, then one row per frame, keeping the
   // selected frames. Reading stops at the end of the file; a row with no
   // volume has a volume of 0.
   std::ifstream ibar(barfile);
   auto read_block = [&](int start, int stop, int step, double& temp, std::string& title, std::vector<double>& u0,
                        std::vector<double>& u1, std::vector<double>& vol) {
      int nrow = 0;
      temp = 0;
      title = "";
      std::string line;
      if (std::getline(ibar, line)) {
         std::istringstream iss(line);
         iss >> nrow >> temp;
         std::getline(iss >> std::ws, title);
         title.erase(title.find_last_not_of(" \t\r") + 1);
      }
      for (int i = 1; i <= nrow and std::getline(ibar, line); ++i) {
         if (i >= start and (i - start) % step == 0 and i <= stop) {
            double e0 = 0, e1 = 0, v = 0;
            std::sscanf(line.c_str(), "%*d %lf %lf %lf", &e0, &e1, &v);
            u0.push_back(e0);
            u1.push_back(e1);
            vol.push_back(v);
         }
      }
   };

   // read potential energies and volumes for both trajectories
   BarData d;
   double tempa, tempb;
   std::string titlea, titleb;
   read_block(starta, stopa, stepa, tempa, titlea, d.ua0, d.ua1, d.vola);
   read_block(startb, stopb, stepb, tempb, titleb, d.ub0, d.ub1, d.volb);
   ibar.close();
   int nfrma = d.vola.size(), nfrmb = d.volb.size();

   // provide info about trajectories and number of frames
   print(out,
      "\n"
      " Simulation Trajectory A and Thermodynamic State 0 :\n"
      "\n"
      " %s\n",
      titlea);
   print(out,
      " Number of Frames :    %8d\n"
      " Temperature :       %10.2lf\n",
      nfrma, tempa);
   print(out,
      "\n"
      " Simulation Trajectory B and Thermodynamic State 1 :\n"
      "\n"
      " %s\n",
      titleb);
   print(out,
      " Number of Frames :    %8d\n"
      " Temperature :       %10.2lf\n",
      nfrmb, tempb);

   // stop if either trajectory has no frames selected
   for (auto t : {std::make_pair('A', nfrma), std::make_pair('B', nfrmb)}) {
      if (t.second == 0) {
         print(out,
            "\n"
            " BAR  --  No Frames Selected from Trajectory %c\n",
            t.first);
         std::fflush(out);
         TINKER_THROW(format("BAR  --  No Frames Selected from Trajectory %c", t.first));
      }
   }

   // stop if a trajectory has no temperature
   for (auto t : {std::make_pair('A', tempa), std::make_pair('B', tempb)}) {
      if (not(t.second > 0)) {
         print(out,
            "\n"
            " BAR  --  No Temperature Found for Trajectory %c\n",
            t.first);
         std::fflush(out);
         TINKER_THROW(format("BAR  --  No Temperature Found for Trajectory %c", t.first));
      }
   }

   // stop if the trajectories differ in temperature
   if (std::fabs(tempa - tempb) > 0.01) {
      print(out,
         "\n"
         " BAR  --  Trajectory Temperatures of%8.2lf and%8.2lf K Differ;\n"
         " FEP and BAR Require a Single Temperature\n",
         tempa, tempb);
      std::fflush(out);
      TINKER_THROW("BAR  --  FEP and BAR Require a Single Temperature");
   }

   // the PV terms need a volume for every frame of both trajectories, or
   // for none; a missing volume would enter the averages as zero
   auto count_vol = [](const std::vector<double>& vol) {
      return (int)std::count_if(vol.begin(), vol.end(), [](double v) { return v != 0; });
   };
   const int nvola = count_vol(d.vola), nvolb = count_vol(d.volb);
   const bool usevol = nvola + nvolb > 0;
   if (usevol) {
      for (auto t : {std::make_tuple('A', nvola, nfrma), std::make_tuple('B', nvolb, nfrmb)}) {
         if (std::get<1>(t) < std::get<2>(t)) {
            print(out,
               "\n"
               " BAR  --  Volumes Missing for%8d of%8d Frames of Trajectory %c;\n"
               " PV Terms Need a Volume for Every Frame\n",
               std::get<2>(t) - std::get<1>(t), std::get<2>(t), std::get<0>(t));
            std::fflush(out);
            TINKER_THROW(format("BAR  --  Volumes Missing for Trajectory %c", std::get<0>(t)));
         }
      }
   }

   // get the pressure for the PV terms if volumes are present
   double patm = 1;
   if (usevol) {
      // only a missing or negative argument asks for the pressure, and any
      // value entered is accepted
      patm = -1;
      nextarg(string, exist);
      if (exist)
         ioReadString(patm, string);
      bool asked = false;
      ioReadStream(patm,
         "\n"
         " Enter Pressure in Atmospheres [1.0] :  ",
         1.0, [&](double p) {
            if (asked)
               return false;
            asked = true;
            return p < 0;
         });
      if (patm <= 0)
         patm = 1;
      print(out,
         "\n"
         " Pressure (Atm) :    %10.2lf\n",
         patm);
   }

   // set the number of bootstrap trials to be generated
   int nfrm = std::max(nfrma, nfrmb);
   int nbst = std::max(2, std::min(100000, (int)std::round(1.e8 / nfrm)));

   // find the free energy, enthalpy and entropy via FEP and BAR
   barPair(d, tempa, tempb, patm, nbst);
}
}

namespace tinker {
void xBar(int, char**)
{
   initial();

   char string[MAX_NCHAR];
   auto out = stdout;
   bool exist;
   int mode = 0;
   auto invalid_mode = [](int m) { return m != 1 and m != 2; };
   const char* mode_string1 = R"(
 The Tinker Thermodynamic Perturbation Utility Can :

    (1) Create BAR File with Perturbed Potential Energies
    (2) Compute Thermodynamic Values from Tinker BAR File
)";

   nextarg(string, exist);
   if (exist)
      ioReadString(mode, string);
   if (invalid_mode(mode))
      print(out, "%s", mode_string1);
   ioReadStream(mode,
      "\n"
      " Enter the Number of the Desired Choice :  ",
      0, invalid_mode);

   if (mode == 1)
      xBarMake();
   if (mode == 2)
      xBarCalc();

   tinker_f_final();
}
}
