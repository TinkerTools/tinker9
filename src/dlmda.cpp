#include "ff/dlmda.h"
#include "ff/ost.h"
#include "ff/eost.h"
#include "ff/eabf.h"
#include "ff/ethrmint.h"
#include "ff/efep.h"
#include "ff/amoeba/induce.h"
#include "ff/atom.h"
#include "ff/egvop.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/evdw.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "math/const.h"
#include "math/random.h"
#include "tool/darray.h"
#include "tool/error.h"
#include "tool/externfunc.h"
#include "tool/ioprint.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/dlmda.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/mutant.hh>
#include <tinker/detail/units.hh>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace tinker {
// The forces of the step, which the gradient of an orthogonal space tempering
// trial of lmdaMonte overwrites; allocated at the first such trial.
static grad_prec *gx_lmdamc, *gy_lmdamc, *gz_lmdamc;

void dlmdaData(RcOp op)
{
   // for whichever system is going away
   if ((op & RcOp::DEALLOC) and gx_lmdamc) {
      darray::deallocate(gx_lmdamc, gy_lmdamc, gz_lmdamc);
      gx_lmdamc = nullptr;
      gy_lmdamc = nullptr;
      gz_lmdamc = nullptr;
   }

   if (!use_rel)
      return;

   if (op & RcOp::DEALLOC)
      darray::deallocate(rdt_group);

   if (op & RcOp::ALLOC)
      darray::allocate(n, &rdt_group);

   if (op & RcOp::INIT) {
      std::vector<int> group(n);
      for (int i = 0; i < n; ++i)
         group[i] = mutant::mutg[i];
      darray::copyin(g::q0, n, rdt_group, group.data());
      waitFor(g::q0);
   }
}

void dlmdaData2(RcOp op)
{
   if (op & RcOp::INIT) {
      bool lambda_dynamics = use_dlmda or use_ost or use_meta or use_ti or use_abf or use_fep or use_epdt //
         or use_rel;
      if (lambda_dynamics and not(pltfm_config & Platform::CUDA))
         TINKER_THROW("LAMBDA  --  Lambda dynamics requires the CUDA platform");
   }
}

Lmdamap lmdamapFrom(const char* s)
{
   // s is a Fortran character*3 buffer (no null terminator).
   if (std::strncmp(s, "EXP", 3) == 0)
      return Lmdamap::EXP;
   if (std::strncmp(s, "INV", 3) == 0)
      return Lmdamap::INV;
   if (std::strncmp(s, "APM", 3) == 0)
      return Lmdamap::APM;
   if (std::strncmp(s, "QNT", 3) == 0)
      return Lmdamap::QNT;
   return Lmdamap::QNT;
}

LmdaThMap lmdaThMapFrom(const char* s)
{
   // s is a Fortran character*3 buffer (no null terminator).
   if (std::strncmp(s, "SIN", 3) == 0)
      return LmdaThMap::SIN;
   return LmdaThMap::TRI;
}

RelStage relStageFrom(const char* s)
{
   // s is a Fortran character*4 buffer (no null terminator).
   if (std::strncmp(s, "LIG1", 4) == 0)
      return RelStage::LIG1;
   if (std::strncmp(s, "LIG2", 4) == 0)
      return RelStage::LIG2;
   return RelStage::VDWM;
}

// dlambda.f:grpscale. The environment is unscaled and the charging ligand
// carries the coupling value; a staged relative leg annihilates the other
// ligand, and the LIG2 leg charges ligand 2 instead of ligand 1.
GrpScale grpScale(double lmda)
{
   GrpScale c;
   c.s[0] = 1;
   c.ds[0] = 0;
   c.s[1] = lmda;
   c.ds[1] = 1;
   c.s[2] = 0;
   c.ds[2] = 0;
   if (use_rel and relstage == RelStage::LIG2) {
      c.s[1] = 0;
      c.ds[1] = 0;
      c.s[2] = lmda;
      c.ds[2] = 1;
   }
   return c;
}

CtScale ctScale()
{
   const bool driven = pdlmdaActive();
   double w = 0, b = 0, d2 = 0;
   GrpScale g0, g1, g;
   if (use_epdt) {
      double dw, d2w;
      dtWeight(plam, epdtexp, w, dw, d2w);
      b = driven ? dw * dpldlmda : 0;
      d2 = driven ? d2w * dpldlmda * dpldlmda + dw * d2pldlmda2 : 0;
      g0 = grpScale(0);
      g1 = grpScale(1);
   } else {
      g = grpScale(plam);
   }

   CtScale c;
   c.live = 0;
   for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
         double sc, dl, d2l;
         if (use_epdt) {
            double s0 = g0.s[i] * g0.s[k];
            double s1 = g1.s[i] * g1.s[k];
            sc = s0 + w * (s1 - s0);
            dl = b * (s1 - s0);
            d2l = d2 * (s1 - s0);
         } else {
            double dsc = g.ds[i] * g.s[k] + g.s[i] * g.ds[k];
            double d2sc = 2 * g.ds[i] * g.ds[k];
            sc = g.s[i] * g.s[k];
            dl = driven ? dsc * dpldlmda : 0;
            d2l = driven ? d2sc * dpldlmda * dpldlmda + dsc * d2pldlmda2 : 0;
         }
         c.sc[i][k] = sc;
         c.dl[i][k] = dl;
         c.d2l[i][k] = d2l;
         if (sc != 0 or dl != 0 or d2l != 0)
            c.live |= 1u << (3 * i + k);
      }
   }
   return c;
}

const int* emGroup()
{
   return use_rel ? rdt_group : mut;
}

RdtMask chargeMask()
{
   return (use_rel and relstage == RelStage::LIG2) ? RdtMask::LIGB : RdtMask::LIGA;
}

RdtMask coupledMask()
{
   return static_cast<RdtMask>(static_cast<unsigned>(RdtMask::ENV) | static_cast<unsigned>(chargeMask()));
}

RdtMask elecLiveMask()
{
   // The charging ligand is dead as well when it is held uncharged with no
   // lambda derivative left to take: the VDWM leg, or an absolute vdw leg.
   bool dead = elam == 0 and plam == 0 and not edlmdaActive() and not pdlmdaActive();
   return dead ? RdtMask::ENV : coupledMask();
}

const int* elecDeadGroups()
{
   // A relative leg always annihilates one ligand; otherwise only a dead
   // charging ligand leaves a group out of the mask.
   bool anydead = use_rel or elecLiveMask() == RdtMask::ENV;
   return anydead ? emGroup() : nullptr;
}

// E = w*E1 + (1-w)*E0, so dE/dl = dw*(E1-E0) and d2E/dl2 = d2w*(E1-E0), which
// is what the endpoint mixing used to form from two finished endpoints. Both terms of
// the second-derivative chain rule are proportional to E1-E0, so the whole of
// it collapses onto one scalar per endpoint. need0/need1 are not consulted: a
// dead endpoint already falls out with every weight at zero.
void dtWeightsToCoef(DtCoef& c, double w, double dw, double d2w, double chain, double d2chain, bool driven)
{
   c.a0 = 1 - w;
   c.a1 = w;
   double b = driven ? dw * chain : 0;
   double d2 = driven ? d2w * chain * chain + dw * d2chain : 0;
   c.b0 = -b;
   c.b1 = b;
   c.c0 = -d2;
   c.c1 = d2;
}

void dtWeight(double x, int nexp, double& w, double& dw, double& d2w)
{
   // The linear case is taken separately so that a zero sub-lambda never
   // reaches a zero power.
   w = std::pow(x, nexp);
   dw = 0.0;
   d2w = 0.0;
   if (nexp == 1) {
      dw = 1.0;
   } else if (nexp >= 2) {
      dw = (double)nexp * std::pow(x, nexp - 1);
      d2w = (double)nexp * (double)(nexp - 1) * std::pow(x, nexp - 2);
   }
}

void dtNeed(double w, double dw, double d2w, double chain, double d2chain, bool& need0, bool& need1)
{
   double c1 = dw * chain;
   double c2 = d2w * chain * chain + dw * d2chain;
   need1 = (w != 0.0) or (c1 != 0.0) or (c2 != 0.0);
   need0 = (w != 1.0) or (c1 != 0.0) or (c2 != 0.0);
}

void dtWeightNeed(double sublmda, int dtexp, double chain, double d2chain, //
   double& w, double& dw, double& d2w, bool& need0, bool& need1)
{
   dtWeight(sublmda, dtexp, w, dw, d2w);
   dtNeed(w, dw, d2w, chain, d2chain, need0, need1);
}

void LmdaBuffer::manage(RcOp op, int flag, bool driven, EnergyBuffer* dl1, EnergyBuffer* dl2,
   energy_prec* term1, energy_prec* term2)
{
   const bool rc_a = flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      if (rc_a and mDl1)
         darray::deallocate(*mDl1, *mDl2);
      if (mDl1) {
         *mDl1 = nullptr;
         *mDl2 = nullptr;
      }
      mDl1 = nullptr;
      mDl2 = nullptr;
      mTerm1 = nullptr;
      mTerm2 = nullptr;
      mFlag = 0;
      mMask = 0;
      mDriven = false;
   }

   if (op & RcOp::ALLOC) {
      mDl1 = dl1;
      mDl2 = dl2;
      mTerm1 = term1;
      mTerm2 = term2;
      mFlag = flag;
      mDriven = driven;

      const int mask = lmdaDerivMask(flag, driven);
      mMask = mask;
      *dl1 = nullptr;
      *dl2 = nullptr;
      if (mask & calc::energy_dlmda1) {
         if (rc_a)
            darray::allocate(bufferSize(), dl1);
         else
            *dl1 = dedl_buf;
      }
      if (mask & calc::energy_dlmda2) {
         if (rc_a)
            darray::allocate(bufferSize(), dl2);
         else
            *dl2 = d2edl2_buf;
      }
   }
}

void LmdaBuffer::zero(int vers) const
{
   if (not(mFlag & calc::analyz) or not mDl1)
      return;
   const int mask = lmdaDerivVers(vers, mDriven) & mMask;
   if (mask & calc::energy_dlmda1)
      darray::zero(g::q0, bufferSize(), *mDl1);
   if (mask & calc::energy_dlmda2)
      darray::zero(g::q0, bufferSize(), *mDl2);
}

void LmdaBuffer::flush(int vers) const
{
   if (not(mFlag & calc::analyz) or not mDl1)
      return;
   const int mask = lmdaDerivVers(vers, mDriven) & mMask;
   if (mask & calc::energy_dlmda1) {
      energy_prec e = energyReduce(*mDl1);
      *mTerm1 += e;
      dedl += e;
   }
   if (mask & calc::energy_dlmda2) {
      energy_prec e = energyReduce(*mDl2);
      *mTerm2 += e;
      d2edl2 += e;
   }
}

static LmdaBuffer ep_dl;

void polarLmdaData(RcOp op)
{
   if (not use(Potent::POLAR) and not use(Potent::CHGTRN))
      return;
   ep_dl.manage(op, rc_flag, use_pdlmda, &depdl_buf, &d2epdl2_buf, &depdl, &d2epdl2);
}

void polarLmdaZero(int vers)
{
   ep_dl.zero(vers);
}

void polarLmdaFlush(int vers)
{
   ep_dl.flush(vers);
}

int dtPassList(DtPass out[2])
{
   // Endpoint 0 zeroes both ligands and endpoint 1 keeps the charging one, with
   // the other ligand of a staged relative leg annihilated (dlambda.f:grpscale);
   // there is nothing the two share.
   out[0] = {RdtMask::ENV, true, false};
   out[1] = {coupledMask(), false, true};
   return 2;
}

void dtPassWeights(const DtCoef& c, const DtPass& p, real& wa, real& wb, real& wc)
{
   wa = (p.in0 ? c.a0 : 0) + (p.in1 ? c.a1 : 0);
   wb = (p.in0 ? c.b0 : 0) + (p.in1 ? c.b1 : 0);
   wc = (p.in0 ? c.c0 : 0) + (p.in1 ? c.c1 : 0);
}


void dlmda_mech()
{
   lambda = mutant::lambda;

   use_dlmda = dlmda::use_dlmda;
   use_d2lmda = dlmda::use_d2lmda;
   use_epdt = dlmda::use_epdt;
   use_prst = dlmda::use_prst;
   use_mainlmda = dlmda::use_mainlmda;
   use_rel = mutant::use_rel;

   use_edlmda = (dlmda::use_edlmda != 0);
   use_pdlmda = (dlmda::use_pdlmda != 0);
   use_vdlmda = (dlmda::use_vdlmda != 0);

   use_epast = use_prst and use_pdlmda;

   // Every multipole run takes the single topology path; a relative one scales
   // each ligand group on its own (grpScale).
   use_emast = use_edlmda;

   epdtexp = dlmda::epdtexp;

   // which lambda-dynamics method owns the main lambda.
   use_ost = dlmda::use_ost;
   use_meta = dlmda::use_meta;
   use_ti = dlmda::use_ti;
   use_fep = dlmda::use_fep;
   use_abf = dlmda::use_abf;
   dlmda::use_ostdyn = use_ost;
   dlmda::use_metadyn = use_meta;
   dlmda::use_abfdyn = use_abf;

   // adaptive lambda bias state shared by OST and ABF.
   lmdastep = dlmda::lmdastep;
   lmdaintv = dlmda::lmdaintv;
   lmdanpa = dlmda::lmdanpa;
   lmdanpb = dlmda::lmdanpb;
   lmdanpc = dlmda::lmdanpc;
   nlmda = dlmda::nlmda;
   nlmdahist = 0;
   sizelmdahist = 0;
   wlmda = dlmda::wlmda;
   wlmda2 = dlmda::wlmda2;
   lmdaparatio = dlmda::lmdaparatio;
   lmdapbratio = dlmda::lmdapbratio;
   lmdatheta = dlmda::lmdatheta;
   lmdavtheta = dlmda::lmdavtheta;
   lmdathmap = lmdaThMapFrom(dlmda::lmdathmap);
   lmdathalpha = dlmda::lmdathalpha;
   lmdamass = dlmda::lmdamass;
   lmdafric = dlmda::lmdafric;
   lmdadt = dlmda::lmdadt;
   use_lmdacv = (dlmda::use_lmdacv != 0);
   lmdacvstd = dlmda::lmdacvstd;
   lmdacvrat = dlmda::lmdacvrat;
   use_lmdamc = (dlmda::use_lmdamc != 0);
   lmdamctrial = dlmda::lmdamctrial;
   lmdamcmove = dlmda::lmdamcmove;
   nlmdamctry = 0;
   nlmdamcacc = 0;

   // dedl is owned by the energy routines and zeroed by zeroEGV each step.
   deffdl = 0;
   lmdaddgdl = 0;
   lmdaavg = 0;
   lmdastd = 0;
   dedlavg = 0;
   dedlstd = 0;
   lmdadeltag = 0;

   elmdaexp = dlmda::elmdaexp;
   plmdaexp = dlmda::plmdaexp;
   vlmdaexp = dlmda::vlmdaexp;
   elmdainvn = dlmda::elmdainvn;
   plmdainvn = dlmda::plmdainvn;
   vlmdainvn = dlmda::vlmdainvn;
   elmdainveps = dlmda::elmdainveps;
   plmdainveps = dlmda::plmdainveps;
   vlmdainveps = dlmda::vlmdainveps;
   elmdaapmn = dlmda::elmdaapmn;
   plmdaapmn = dlmda::plmdaapmn;
   vlmdaapmn = dlmda::vlmdaapmn;
   elmdaapmrho = dlmda::elmdaapmrho;
   plmdaapmrho = dlmda::plmdaapmrho;
   vlmdaapmrho = dlmda::vlmdaapmrho;

   elmdamap = lmdamapFrom(dlmda::elmdamap);
   plmdamap = lmdamapFrom(dlmda::plmdamap);
   vlmdamap = lmdamapFrom(dlmda::vlmdamap);

   use_elmdamap = (dlmda::use_elmdamap != 0);
   use_plmdamap = (dlmda::use_plmdamap != 0);
   use_vlmdamap = (dlmda::use_vlmdamap != 0);

   deldlmda = dlmda::deldlmda;
   dpldlmda = dlmda::dpldlmda;
   dvldlmda = dlmda::dvldlmda;
   d2eldlmda2 = 0.0;
   d2pldlmda2 = 0.0;
   d2vldlmda2 = 0.0;

   qntelmda0 = dlmda::qntelmda0;
   qntelmda1 = dlmda::qntelmda1;
   qntplmda0 = dlmda::qntplmda0;
   qntplmda1 = dlmda::qntplmda1;
   qntvlmda0 = dlmda::qntvlmda0;
   qntvlmda1 = dlmda::qntvlmda1;

   relstage = relStageFrom(dlmda::relstage);
   for (int i = 0; i < 3; ++i) {
      use_relstg[i] = (dlmda::use_relstg[i] != 0);
      relstglo[i] = dlmda::relstglo[i];
      relstghi[i] = dlmda::relstghi[i];
   }
}

void avgstd(const std::vector<double>& v, int begin, int count, double& avg, double& sd)
{
   if (count < 1) {
      avg = 0.0;
      sd = 0.0;
      return;
   }

   // Shifted mean, to keep the sum of squares well conditioned.
   const double k = v[begin];
   double total = 0.0, totalsq = 0.0;
   for (int i = begin; i < begin + count; ++i) {
      double d = v[i] - k;
      total += d;
      totalsq += d * d;
   }
   avg = k + total / (double)count;
   double var = (totalsq - total * total / (double)count) / (double)count;
   sd = std::sqrt(var > 0.0 ? var : 0.0);
}

// depcriteria -- whether the samples of the last interval are converged enough
// to be added to the lambda bias, by comparing their deviation against a
// tolerance with both absolute and relative parts; every interval is accepted
// unless the convergence gate is turned on (dlambda.f:depcriteria).
bool depcriteria(double avg, double sd)
{
   // accept every interval when the gate is off
   if (not use_lmdacv)
      return true;

   // accept only a deviation strictly inside the tolerance
   double tolerance = lmdacvstd + lmdacvrat * std::abs(avg);
   return tolerance > 0.0 && sd / tolerance < 1.0;
}

// Lambda bin index of a lambda value, clamped to the grid (dlambda.f:lmdabin).
int lmdaBin(double lambda)
{
   int b = (int)std::lround(lambda / wlmda) + 1;
   if (b < 1)
      b = 1;
   if (b > nlmda)
      b = nlmda;
   return b;
}

// Divides the lambda sample interval into the phase that propagates the lambda
// particle, the phase that equilibrates at the frozen lambda and the phase that
// averages dU/dlambda at that same fixed lambda, the averaging phase taking
// whatever the other two leave (dlambda.f:setlmdaphase).
void setLmdaPhase()
{
   // divide the interval, keeping at least one propagation step and at least
   // two samples to average
   if (lmdaintv < 1)
      lmdaintv = 1;
   lmdanpa = (int)(lmdaparatio * (double)lmdaintv);
   lmdanpb = (int)(lmdapbratio * (double)lmdaintv);
   lmdanpa = std::max(1, std::min(lmdanpa, lmdaintv - 1));
   lmdanpb = std::max(0, std::min(lmdanpb, lmdaintv - lmdanpa));
   lmdanpc = lmdaintv - lmdanpa - lmdanpb;
   while (lmdanpc < 2 and lmdanpb > 0) {
      lmdanpb -= 1;
      lmdanpc += 1;
   }
   while (lmdanpc < 2 and lmdanpa > 1) {
      lmdanpa -= 1;
      lmdanpc += 1;
   }
}

// Maps theta onto the main lambda, as sine squared or a smoothed triangle, and
// returns dlambda/dtheta. The triangle map is
//    lambda = 1/2 - asin(a*cos(2*theta)) / (2*asin(a))
// which becomes sin(theta)^2 as the sharpness a goes to zero and a triangle
// wave, flat in lambda, as a goes to one. Both maps have period pi in theta
// (dlambda.f:lmdathetamap).
void lmdaThetaMap(double theta, double& lmda, double& dldth)
{
   if (lmdathmap == LmdaThMap::TRI) {
      double alpha = lmdathalpha;
      double asina = std::asin(alpha);
      double costh = std::cos(2.0 * theta);
      lmda = 0.5 - std::asin(alpha * costh) / (2.0 * asina);
      dldth = alpha * std::sin(2.0 * theta) / (asina * std::sqrt(1.0 - alpha * alpha * costh * costh));
   } else {
      double sinth = std::sin(theta);
      lmda = sinth * sinth;
      dldth = std::sin(2.0 * theta);
   }
}

// Inverts the theta map onto the principal branch [0, pi/2]
// (dlambda.f:lmdathetainv).
void lmdaThetaInv(double lmda, double& theta)
{
   double lclip = std::min(1.0, std::max(0.0, lmda));
   if (lmdathmap == LmdaThMap::TRI) {
      double alpha = lmdathalpha;
      double asina = std::asin(alpha);
      double arg = std::sin((1.0 - 2.0 * lclip) * asina) / alpha;
      arg = std::min(1.0, std::max(-1.0, arg));
      theta = 0.5 * std::acos(arg);
   } else {
      theta = std::asin(std::sqrt(lclip));
   }
}

// BAOAB Langevin propagation of the auxiliary lambda particle in theta space,
// where lambda follows the theta map lmdathmap; it draws the random kick and
// has lmdaLangStep take the step (dlambda.f:lmdalangevin).
void lmdaLangevin()
{
   if (lmdadt <= 0.0)
      return;
   if (lmdamass <= 0.0)
      return;

   // draw the random kick only when the friction puts it to use
   double gauss = 0.0;
   if (lmdafric > 0.0)
      gauss = normal<double>();
   lmdaLangStep(gauss);
}

// One Langevin step of the lambda particle in theta space, with the standard
// normal deviate gauss as its random kick. The force on theta is in kcal/mole,
// so the mass is in kcal/mole*ps**2 and the size of the kick follows from the
// thermal energy in kcal/mole (dlambda.f:lmdalangstep).
void lmdaLangStep(double gauss)
{
   // force on theta from dU/dlambda and the theta map derivative
   double lmda, dldth;
   lmdaThetaMap(lmdatheta, lmda, dldth);
   double force = -deffdl * dldth;
   double gamma = std::max(0.0, lmdafric);
   if (gamma > 0.0) {
      double c = std::exp(-gamma * lmdadt);
      double ktm = units::gasconst * bath::kelvin / lmdamass;
      double sigma = std::sqrt(ktm * (1.0 - c * c));
      lmdavtheta = c * lmdavtheta + (1.0 - c) * force / (gamma * lmdamass) + sigma * gauss;
   } else {
      lmdavtheta = lmdavtheta + lmdadt * force / lmdamass;
   }

   // update theta and wrap it into the periodic interval [0,pi). The period is
   // M_PI and not tinker::pi, which is only a float in a mixed precision build.
   constexpr double dpi = M_PI;
   lmdatheta = lmdatheta + lmdadt * lmdavtheta;
   lmdatheta = lmdatheta - dpi * std::floor(lmdatheta / dpi);
   if (lmdatheta >= dpi)
      lmdatheta = 0.0;
   if (lmdatheta < 0.0)
      lmdatheta = 0.0;

   // map theta back to the main lambda
   lmdaThetaMap(lmdatheta, lambda, dldth);
}

// elmdaDyn -- advances the adaptive lambda bias by one dynamics step for
// orthogonal space tempering, metadynamics or adaptive biasing force; it builds
// the effective lambda derivative, saves the interval samples, deposits the
// interval average at the end of each interval and propagates the lambda
// particle; when the lambda is moved by Monte Carlo trials the particle is left
// alone and lmdaMonte changes the lambda once the step is complete
// (dlambda.f:elmdadyn).
void elmdaDyn(int istep)
{
   int im = istep % lmdaintv;
   int isamp = (istep - 1) % lmdaintv;

   // build the effective lambda derivative from the unbiased dedl left behind
   // by the energy call and the bias saved this step by eostBias or eabfBias
   if (use_ost) {
      ostdgdl = bdgdl + bdgdfl * d2edl2;
      lmdaddgdl = lmdadfdl;
      deffdl = dedl + ostdgdl - lmdaddgdl;
   } else if (use_meta) {
      deffdl = dedl + bdgdl;
   } else if (use_abf) {
      lmdaddgdl = lmdadfdl;
      deffdl = dedl - lmdaddgdl;
   }

   // save all values in the interval, but average only after the propagation
   // and equilibration phases; only OST and ABF keep dU/dlambda, which they
   // average at a fixed lambda
   lmdallist[isamp] = lambda;
   if (use_ost or use_abf)
      lmdaflist[isamp] = dedl;

   // deposit the interval average every lmdaintv steps; an OST or ABF interval
   // is kept only when its samples are converged
   if (im == 0) {
      int nskip = lmdanpa + lmdanpb;
      avgstd(lmdallist, nskip, lmdanpc, lmdaavg, lmdastd);
      if (use_ost or use_abf) {
         avgstd(lmdaflist, nskip, lmdanpc, dedlavg, dedlstd);
         if (depcriteria(dedlavg, dedlstd)) {
            if (use_ost)
               ostDeposit(istep);
            if (use_abf)
               abfDeposit(istep);
         }
      } else if (use_meta) {
         metaDeposit(istep);
      }
   }

   // Monte Carlo trials move the lambda in place of the lambda particle
   if (use_lmdamc)
      return;

   // propagate the lambda particle; OST and ABF hold lambda fixed after the
   // propagation phase
   if (use_ost or use_abf) {
      if (isamp < lmdanpa)
         lmdaLangevin();
   } else {
      lmdaLangevin();
   }
}

// Free energy at the current lambda by piecewise-linear integration of the mean
// force of the lambda bins, and its lambda derivative (dlambda.f:efreelmda).
void efreeLmda(double& eflmda, double& dfdl)
{
   eflmda = 0;
   dfdl = 0;
   if (lambda <= 0.0) {
      dfdl = lmdafmean[1];
      return;
   }
   for (int il0 = 1; il0 <= nlmda - 1; ++il0) {
      int il1 = il0 + 1;
      double lmda0 = (double)(il0 - 1) * wlmda;
      double lmda1 = (double)(il1 - 1) * wlmda;
      double fl0 = lmdafmean[il0];
      double fl1 = lmdafmean[il1];
      double slope = (fl1 - fl0) / wlmda;
      if (lambda <= lmda1) {
         double xx = lambda - lmda0;
         eflmda += fl0 * xx + 0.5 * slope * xx * xx;
         dfdl = fl0 + slope * xx;
         return;
      }
      eflmda += 0.5 * (fl0 + fl1) * wlmda;
   }
   dfdl = lmdafmean[nlmda];
}

// Total free energy change by trapezoid integration of the mean force of the
// lambda bins (dlambda.f:efreetot).
double efreeTot()
{
   double tot = 0;
   for (int il = 1; il <= nlmda - 1; ++il)
      tot += 0.5 * (lmdafmean[il] + lmdafmean[il + 1]) * wlmda;
   return tot;
}

void lmdaWinData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
      lmdawinend.clear();
      lmdawinfrac.clear();
      lmdawinlist.clear();
   }

   if (op & RcOp::INIT) {
      lmdawin = dlmda::lmdawin;
      lmdawineq = dlmda::lmdawineq;
      lmdawinlen = dlmda::lmdawinlen;
      nlmdawin = dlmda::nlmdawin;
      lmdawinratio = dlmda::lmdawinratio;

      // The boundaries wait for the length of the run; lmdaWinFromFortran
      // reads them once the windows are laid out.
      lmdawinend.clear();

      // Only a sequential lambda method has a schedule to adopt.
      if (use_ti or use_fep) {
         lmdawinfrac.assign(dlmda::lmdawinfrac, dlmda::lmdawinfrac + nlmdawin);
         lmdawinlist.assign(dlmda::lmdawinlist, dlmda::lmdawinlist + nlmdawin);
      } else {
         lmdawinfrac.clear();
         lmdawinlist.clear();
      }
   }
}

void lmdaWinFromFortran()
{
   nlmdawin = dlmda::nlmdawin;
   lmdawin = dlmda::lmdawin;
   lmdawineq = dlmda::lmdawineq;
   lmdawinlen = dlmda::lmdawinlen;
   lmdawinend.assign(dlmda::lmdawinend, dlmda::lmdawinend + nlmdawin);
   lambda = mutant::lambda;
}

void setLmdaWin()
{
   lmdaWinSize(lmdawin, lmdawinlen, lmdawineq);
}

void lmdaWinSize(int iwin, int& nlen, int& neq)
{
   // Measure the window against the preceding boundary.
   nlen = lmdawinend[iwin - 1];
   if (iwin > 1)
      nlen = lmdawinend[iwin - 1] - lmdawinend[iwin - 2];
   neq = (int)((double)nlen * lmdawinratio);
}

void lmdaWinPhase(int istep, int& iprod)
{
   // No step is in production outside of the window schedule.
   iprod = 0;
   if (lmdawinend.empty())
      return;
   if (lmdawin < 1 or lmdawin > nlmdawin)
      return;

   // Count the steps taken since the window finished equilibrating.
   int istart = (lmdawin > 1) ? lmdawinend[lmdawin - 2] : 0;
   iprod = std::max(0, istep - istart - lmdawineq);
}

// The bias energy at the current lambda of the methods whose bias follows
// lambda alone, adaptive biasing force and metadynamics, from the history as it
// now stands. No saved value is changed (dlambda.f:elmdabias).
static double elmdaBias()
{
   double vbias, slope;
   if (use_abf) {
      efreeLmda(vbias, slope);
      return -vbias;
   }
   eMetaBias(lambda, vbias, slope);
   return vbias;
}

void lmdaMonte(int istep)
{
   // OST and ABF move lambda only in the propagation phase of the sample
   // interval, while metadynamics moves it at any step
   const int isamp = (istep - 1) % lmdaintv;
   if ((use_ost or use_abf) and isamp >= lmdanpa)
      return;

   // decide whether to attempt a lambda change at this step
   if (random<double>() >= 1.0 / lmdamctrial)
      return;

   // draw the trial lambda and reflect it into the lambda range
   const double kt = units::gasconst * bath::kelvin;
   double ltrial = lambda + lmdamcmove * (2 * random<double>() - 1);
   if (ltrial < 0)
      ltrial = -ltrial;
   if (ltrial > 1)
      ltrial = 2 - ltrial;
   ++nlmdamctry;

   // The energy of the step no longer holds once the coordinates or the box
   // have been set since it was found, as they are by a scaling barostat; nor
   // does the OST energy hold once this step has added to the history.
   const bool moved = esum_epoch != xyz_epoch;
   const bool stale = moved or (use_ost and istep % lmdaintv == 0);

   // Save the native state before any energy is found, so that a rejected move
   // restores just what the step has left. The trial solves start from the
   // direct guess and leave the predictor history to the dynamics.
   const double lold = lambda;
   const energy_prec estep = esum;
   const auto epoch = esum_epoch;
   const LmdaStepState step;
   InducedTrial trial(InducedTrial::History::BYPASS);
   if (use_ost) {
      if (not gx_lmdamc)
         darray::allocate(n, &gx_lmdamc, &gy_lmdamc, &gz_lmdamc);
      darray::copy(g::q0, n, gx_lmdamc, gx);
      darray::copy(g::q0, n, gy_lmdamc, gy);
      darray::copy(g::q0, n, gz_lmdamc, gz);
   }

   // Get the old energy at the native lambda, using the energy of the step
   // unless it is stale. A bias that follows lambda alone is taken from the
   // history as it now stands, which may hold a deposit made after the forces
   // of the step were found; an energy call does the same.
   const int vers = use_ost ? calc::v4 : calc::energy;
   double eold;
   if (stale) {
      // a rescaled structure has left the neighbor lists behind
      nblistRefresh();
      energy(vers);
      eold = esum;
   } else if (use_ost) {
      eold = estep;
   } else {
      eold = estep - (use_abf ? -bostlmda : bgbias) + elmdaBias();
   }

   // Get the new energy at the trial lambda; the OST bias follows dU/dlambda,
   // so it is found along with the lambda derivative.
   lambda = ltrial;
   energy(vers);
   const double enew = esum;

   // acceptance from the change in potential energy plus bias
   const double de = enew - eold;
   bool accept = de <= 0;
   if (not accept and kt > 0)
      accept = random<double>() < std::exp(-de / kt);

   if (accept) {
      // accept the step, and keep the lambda particle at the new lambda
      ++nlmdamcacc;
      if (use_ost)
         ostdgdl = bdgdl + bdgdfl * d2edl2;
      if (use_ost or use_abf)
         lmdaddgdl = lmdadfdl;
      lmdaThetaInv(lambda, lmdatheta);
      lmdavtheta = 0;
      trial.accept();
   } else {
      // Reject the step, and restore values prior to the trial change. The
      // lambda dependent parameters are scaled again by whatever reads them
      // next, the following energy call or the saved output.
      lambda = lold;
      if (use_mainlmda)
         mapSubLambda();
      step.restore();
      esum = estep;
      // The energy calls of the trial marked the energy as current, which the
      // energy of a step that the system had already left behind is not.
      if (moved)
         esum_epoch = epoch;
   }

   // The integrator keeps the forces of the step for its next half kick.
   if (use_ost) {
      darray::copy(g::q0, n, gx, gx_lmdamc);
      darray::copy(g::q0, n, gy, gy_lmdamc);
      darray::copy(g::q0, n, gz, gz_lmdamc);
   }
}

void lmdaMonteCheck()
{
   const bool bias = dlmda::use_ost or dlmda::use_meta or dlmda::use_abf;
   if (dlmda::use_lmdamc and bias and not bath::isothermal)
      TINKER_THROW("DYNAMIC  --  Monte Carlo Lambda Moves Require a Constant Temperature Ensemble");
}

void lmdaWinDyn(int istep)
{
   // advance the adaptive lambda bias and the lambda particle, now that the
   // step and any Monte Carlo box change are complete
   if (use_ost or use_meta or use_abf)
      elmdaDyn(istep);

   // collect dU/dlambda for thermodynamic integration
   if (use_ti)
      etidyn(istep);

   // find trial lambda energies for free energy perturbation
   if (use_fep)
      efepDyn(istep);

   // try a Monte Carlo move of an adaptively biased lambda
   if (use_lmdamc and (use_ost or use_meta or use_abf))
      lmdaMonte(istep);
}

// The gradients held by the integrator are not updated for the new lambda, so
// the first half step of the velocities in a new window still uses forces from
// the previous window (dlambda.f:nextlmdawin).
void nextLmdaWin()
{
   if (lmdawin > nlmdawin)
      return;
   ++lmdawin;
   if (lmdawin <= nlmdawin) {
      lambda = lmdawinlist[lmdawin - 1];
      setLmdaWin();
   }
}

void lmdaWinStep(int istep)
{
   // Nothing is left to advance outside of the window schedule.
   if (lmdawinend.empty())
      return;
   if (lmdawin < 1 or lmdawin > nlmdawin)
      return;

   // Move on to the next lambda window at the window boundary.
   if (istep == lmdawinend[lmdawin - 1])
      nextLmdaWin();
}

}

namespace tinker {
// Power-law sub-lambda map, lmda = x^exponent.
static void sublmdaExp(double x, int exponent, double& lmda, double& dlmda, double& d2lmda)
{
   lmda = std::pow(x, exponent);
   if (exponent == 1) {
      dlmda = 1.0;
      d2lmda = 0.0;
   } else if (exponent == 2) {
      dlmda = 2.0 * x;
      d2lmda = 2.0;
   } else {
      double expnt = (double)exponent;
      dlmda = expnt * std::pow(x, exponent - 1);
      d2lmda = expnt * (expnt - 1.0) * std::pow(x, exponent - 2);
   }
}

// Shifted inverse-power sub-lambda map onto [0,1] (dlambda.f:sublmdainvpower).
static void sublmdaInvPower(double x, int nn, double eps, double& lmda, double& dlmda, double& d2lmda)
{
   double xval = x;
   if (nn <= 1) {
      lmda = xval;
      dlmda = 1.0;
      d2lmda = 0.0;
      return;
   }
   double shift = eps;
   if (shift <= 0.0)
      shift = 0.1;
   double power = 1.0 / (double)nn;
   double root0 = std::pow(shift, power);
   double denom = std::pow(1.0 + shift, power) - root0;
   double base = xval + shift;
   lmda = (std::pow(base, power) - root0) / denom;
   dlmda = power * std::pow(base, power - 1.0) / denom;
   d2lmda = power * (power - 1.0) * std::pow(base, power - 2.0) / denom;
}

// Normalized asymmetric power sub-lambda map onto [0,1] (dlambda.f:sublmdaapm).
// The map runs from zero to one with slope "rho" at the decoupled end and slope
// one at the coupled end, so "rho" sets how much faster the sub-lambda leaves
// the decoupled end than it arrives at the coupled end, and "n" sets how quickly
// that head start decays.
static void sublmdaApm(double x, int nn, double rho, double& lmda, double& dlmda, double& d2lmda)
{
   double xval = x;
   // A flat slope ratio or a degenerate power gives the identity map, the upper
   // bound holding the normalization off of its pole.
   if (nn < 2 or rho <= 1.0 or rho >= (double)nn) {
      lmda = xval;
      dlmda = 1.0;
      d2lmda = 0.0;
      return;
   }

   double rn = (double)nn;
   double left = rn * (rho - 1.0) / (rn - rho);
   double right = left / rn;
   double denom = 1.0 + right;
   double omx = 1.0 - xval;
   double head = left * (1.0 - std::pow(omx, nn + 1)) / (rn + 1.0);
   double tail = right * std::pow(xval, nn + 1) / (rn + 1.0);
   lmda = (xval + head + tail) / denom;
   dlmda = (1.0 + left * std::pow(omx, nn) + right * std::pow(xval, nn)) / denom;
   d2lmda = rn * (right * std::pow(xval, nn - 1) - left * std::pow(omx, nn - 1)) / denom;
}

// Quintic switching polynomial
void quinticTaper(double x, double cut, double off, double& taper, double& dtaper, double& d2taper)
{
   dtaper = 0.0;
   d2taper = 0.0;
   if (x <= cut) {
      taper = 1.0;
      return;
   } else if (x >= off) {
      taper = 0.0;
      return;
   }

   double rinv = 1.0 / (off - cut);
   double u = (x - cut) * rinv;
   double u2 = u * u;
   double v = 1.0 - u;

   taper = 1.0 - u2 * u * (10.0 - 15.0 * u + 6.0 * u2);
   dtaper = -30.0 * u2 * v * v * rinv;
   d2taper = -60.0 * u * v * (1.0 - 2.0 * u) * rinv * rinv;
}

// Maps one sub-lambda from main lambda to its mapping type.
static void mapOne(double lmda, Lmdamap map, double qnt0, double qnt1, int expExp, int invN, double invEps, int apmN,
   double apmRho, double& value, double& dvalue, double& d2value)
{
   if (map == Lmdamap::EXP) {
      sublmdaExp(lmda, expExp, value, dvalue, d2value);
   } else if (map == Lmdamap::INV) {
      sublmdaInvPower(lmda, invN, invEps, value, dvalue, d2value);
   } else if (map == Lmdamap::APM) {
      sublmdaApm(lmda, apmN, apmRho, value, dvalue, d2value);
   } else { // Lmdamap::QNT
      // quantized map: sublambda = 1 - taper.
      double taper, dtaper, d2taper;
      quinticTaper(lmda, qnt0, qnt1, taper, dtaper, d2taper);
      value = 1.0 - taper;
      dvalue = -dtaper;
      d2value = -d2taper;
   }
}

// Finds the declared leg of the staged relative schedule that holds the main
// lambda and makes it the active leg, returning the local coordinate s running
// from zero to one across that leg and its derivative ds with respect to the
// main lambda; a main lambda outside every declared range holds the end state
// of the adjacent leg, where the derivative is zero (dlambda.f:relstgloc).
static void relStgLoc(double lmda, double& s, double& ds)
{
   // Take the last declared leg starting at or below the main lambda, or the
   // first declared leg when the main lambda lies below all.
   int k = -1;
   for (int i = 0; i < 3; ++i)
      if (use_relstg[i] and (k < 0 or lmda >= relstglo[i]))
         k = i;

   // With no leg declared the active leg is left as it was set, and spans the
   // whole main lambda.
   s = lmda;
   ds = 1.0;
   if (k < 0)
      return;

   relstage = static_cast<RelStage>(k);
   const double lo = relstglo[k];
   const double hi = relstghi[k];
   ds = 1.0 / (hi - lo);
   s = std::min(1.0, std::max(0.0, (lmda - lo) * ds));
   if (lmda >= hi)
      s = 1.0;
   if (lmda < lo or lmda > hi)
      ds = 0.0;
}

// Maps the main lambda onto the sub-lambdas of the staged relative schedule,
// where each declared leg spans its own range of the main lambda
// (dlambda.f:maprelstage):
//
//    LIG2   charge ligand 2 with ligand 1 annihilated, its weight
//             falling as the main lambda rises
//    VDWM   both ligands electrostatically annihilated while van der
//             Waals morphs from ligand 2 onto ligand 1
//    LIG1   charge ligand 1 with ligand 2 annihilated, its weight
//             rising with the main lambda
//
// The leg holding the main lambda is walked on its local coordinate, the term
// it drives following its own map across the whole leg, and polarization
// always carries the electrostatic weight.
static void mapRelStage(double lmda)
{
   double s, ds;
   relStgLoc(lmda, s, ds);

   double eval, vval;

   if (relstage == RelStage::VDWM) {
      // The middle leg holds both ligands uncharged, so electrostatics sits at
      // the reference state and leaves the chain rule while van der Waals
      // morphs across its map.
      eval = 0.0;
      deldlmda = 0.0;
      d2eldlmda2 = 0.0;
      mapOne(s, vlmdamap, 0.0, 1.0, vlmdaexp, vlmdainvn, vlmdainveps, //
         vlmdaapmn, vlmdaapmrho, vval, dvldlmda, d2vldlmda2);
      vlam = vval;
   } else if (relstage == RelStage::LIG1) {
      // The ligand 1 leg charges ligand 1 with ligand 2 annihilated and van
      // der Waals already morphed onto ligand 1.
      mapOne(s, elmdamap, 0.0, 1.0, elmdaexp, elmdainvn, elmdainveps, //
         elmdaapmn, elmdaapmrho, eval, deldlmda, d2eldlmda2);
      vlam = 1.0;
      dvldlmda = 0.0;
      d2vldlmda2 = 0.0;
   } else {
      // The ligand 2 leg discharges ligand 2 as the main lambda rises, so its
      // weight is the complement of the map, with van der Waals still on it.
      mapOne(s, elmdamap, 0.0, 1.0, elmdaexp, elmdainvn, elmdainveps, //
         elmdaapmn, elmdaapmrho, eval, deldlmda, d2eldlmda2);
      eval = 1.0 - eval;
      deldlmda = -deldlmda;
      d2eldlmda2 = -d2eldlmda2;
      vlam = 0.0;
      dvldlmda = 0.0;
      d2vldlmda2 = 0.0;
   }

   // Numerical guard on the map complement.
   elam = std::min(1.0, std::max(0.0, eval));

   // Carry the chain rule from the local coordinate of the leg back to the
   // main lambda.
   deldlmda *= ds;
   d2eldlmda2 *= ds * ds;
   dvldlmda *= ds;
   d2vldlmda2 *= ds * ds;

   // Polarization stages with the multipoles, same state same weight.
   plam = elam;
   dpldlmda = deldlmda;
   d2pldlmda2 = d2eldlmda2;
}

bool lmdaSameValue(double a, double b)
{
   constexpr double eps = 1.0e-6;
   return std::fabs(a - b) <= eps;
}

bool polTracksEle()
{
   // On a staged leg polarization carries the electrostatic weight, whatever
   // map and window it names (dlambda.f:maprelstage).
   if (use_rel)
      return true;
   // One sub-lambda driven and the other frozen: they part company as soon as
   // the main lambda moves off the value they happen to share now.
   if (use_elmdamap != use_plmdamap)
      return false;
   // Neither is driven, so today's values are the only values.
   if (not(use_mainlmda and use_elmdamap))
      return lmdaSameValue(elam, plam);

   if (plmdamap != elmdamap)
      return false;

   if (elmdamap == Lmdamap::EXP) {
      return plmdaexp == elmdaexp;
   } else if (elmdamap == Lmdamap::INV) {
      return plmdainvn == elmdainvn and lmdaSameValue(plmdainveps, elmdainveps);
   } else if (elmdamap == Lmdamap::APM) {
      return plmdaapmn == elmdaapmn and lmdaSameValue(plmdaapmrho, elmdaapmrho);
   } else {
      return lmdaSameValue(qntplmda0, qntelmda0) and lmdaSameValue(qntplmda1, qntelmda1);
   }
}

bool usePoleorig()
{
   bool elec = use(Potent::MPOLE) or use(Potent::POLAR) or use(Potent::REPULS);
   // Free energy perturbation moves the main lambda with no lambda derivative
   // in use, so the multipoles are rescaled from the originals there too.
   return elec and (use_dlmda or use_epdt or use_prst or use_fep);
}

void mapSubLambda()
{
   if (use_rel) {
      mapRelStage(lambda);
      return;
   }

   // Map the main lambda onto each sub-lambda that follows it.
   if (use_plmdamap) {
      double pval;
      mapOne(lambda, plmdamap, qntplmda0, qntplmda1, plmdaexp, plmdainvn, plmdainveps, //
         plmdaapmn, plmdaapmrho, pval, dpldlmda, d2pldlmda2);
      plam = pval;
   }
   if (use_elmdamap) {
      double eval;
      mapOne(lambda, elmdamap, qntelmda0, qntelmda1, elmdaexp, elmdainvn, elmdainveps, //
         elmdaapmn, elmdaapmrho, eval, deldlmda, d2eldlmda2);
      elam = eval;
   }
   if (use_vlmdamap) {
      double vval;
      mapOne(lambda, vlmdamap, qntvlmda0, qntvlmda1, vlmdaexp, vlmdainvn, vlmdainveps, //
         vlmdaapmn, vlmdaapmrho, vval, dvldlmda, d2vldlmda2);
      vlam = vval;
   }

}

}
