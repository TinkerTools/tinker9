#include "ff/amoeba/empole.h"
#include "ff/amoeba/mpolestate.h"
#include "ff/atom.h"
#include "ff/dlmda.h"
#include "ff/elec.h"
#include "ff/hippo/cflux.h"
#include "ff/hippo/erepel.h"
#include "ff/modamoeba.h"
#include "ff/modhippo.h"
#include "ff/pme.h"
#include "ff/potent.h"
#include "tool/darray.h"
#include "tool/externfunc.h"
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/bndstr.hh>
#include <tinker/detail/dlmda.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/mpole.hh>
#include <tinker/detail/mutant.hh>

#include <cassert>
#include <vector>

namespace tinker {
TINKER_FVOID2(acc1, cu1, torque, int, grad_prec*, grad_prec*, grad_prec*);
void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz)
{
   TINKER_FCALL2(acc1, cu1, torque, vers, dx, dy, dz);
}

TINKER_FVOID2(acc1, cu1, torque, int, grad_prec*, grad_prec*, grad_prec*, const real*, const real*,
   const real*, VirialBuffer);
void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz, const real* tqx, const real* tqy,
   const real* tqz, VirialBuffer vbuf)
{
   TINKER_FCALL2(acc1, cu1, torque, vers, dx, dy, dz, tqx, tqy, tqz, vbuf);
}

// What the shared multipole buffers were last built from. A helper below skips
// its kernel when the inputs it would read are the ones it read last time.
namespace {
constexpr unsigned long long NEVER = ~0ull;

// Serial of the contents of pole, advanced by every write to it.
unsigned long long pole_serial = 1;
// The lambda mpoleScale() last scaled pole to, the ligand it charged there,
// and the contents it produced. chkpole inverts pole and poleorig together, so
// it keeps the relation.
double pole_lmda = 0;
RdtMask pole_lmda_mask = RdtMask::LIGA;
// The lambda and the ligand the charge penetration and charge flux parameters
// were last scaled to (mpoleScaleChgpenCflux); unset while they are as the
// setup left them. Charge flux rewrites pole at every step, but never these.
bool param_lmda_set = false;
double param_lmda = 0;
RdtMask param_lmda_mask = RdtMask::LIGA;
unsigned long long pole_lmda_serial = 0;

// The epoch and contents of the last charge flux applied to pole; chkpole
// inverts the y components only, so it keeps the monopoles the flux wrote.
unsigned long long flux_epoch = NEVER, flux_serial = 0;
// The contents of pole the flux was first applied to. The flux rewrites the
// monopoles alone, from their scaled values, so applying it again keeps this.
unsigned long long flux_base_serial = 0;

// Epoch at which chkpole last ran; its result is recorded in zaxis.
unsigned long long chk_epoch = NEVER;

enum class RpoleSrc
{
   NONE,   // nothing recorded; no request carries it, so it never matches
   POLE,   // rotpole(false), from pole
   ORIG,   // rotpole(true), from poleorig
   MASKED, // rotpoleState(), from poleorig
};

struct RpoleKey
{
   RpoleSrc src = RpoleSrc::NONE;
   RdtMask mask = RdtMask::ALL;
   const int* group = nullptr;
   unsigned long long epoch = NEVER, pole = 0;

   bool operator==(const RpoleKey& k) const
   {
      return src == k.src and mask == k.mask and group == k.group and epoch == k.epoch and pole == k.pole;
   }
};
RpoleKey rpole_key;
unsigned long long rpole_serial = 1;

// cmp is a function of rpole alone, or with dlmda also of the lambda scale;
// dl adds dlcmp without changing cmp.
struct CmpKey
{
   bool valid = false, dlmda = false, dl = false;
   double lmda = 0;
   const int* group = nullptr;
   unsigned long long rpole = 0;
};
CmpKey cmp_key;
// Serial of the contents of cmp, advanced by every write to it.
unsigned long long cmp_serial = 1;

// fmp, fphi and cphi hold the reciprocal space potential of the cmp with this
// serial, gathered at this epoch on a PME grid with these parameters. With vir,
// vir_m holds the virial of the same convolution: every convolution into vir_m
// zeroes it first and produces fphi, so the two are always of one convolution.
struct FphiKey
{
   bool valid = false, vir = false;
   unsigned long long cmp = 0, epoch = NEVER;
   PME::Params params{0, 0, 0, 0, 0};
};
FphiKey fphi_key;

// Whether the term that last called mpoleBegin() computes the virial, so that
// the polarization direct field also needs the convolution virial in vir_m.
bool recip_vir = false;

struct SplineKey
{
   unsigned long long epoch = NEVER;
   int level = 0;
};
std::vector<SplineKey> spline_key;
}

void mpoleStateReset()
{
   ++pole_serial;
   param_lmda_set = false;
   chk_epoch = NEVER;
   rpole_key = RpoleKey();
   cmp_key = CmpKey();
   ++cmp_serial;
   fphi_key = FphiKey();
   spline_key.clear();
}

void mpoleCmpClobbered()
{
   cmp_key.valid = false;
   ++cmp_serial;
   fphi_key.valid = false;
}

void mpoleFluxApplied()
{
   if (pole_serial != flux_serial)
      flux_base_serial = pole_serial;
   flux_epoch = xyz_epoch;
   flux_serial = ++pole_serial;
}

bool mpoleFluxCurrent()
{
   return flux_epoch == xyz_epoch and flux_serial == pole_serial;
}

void mpoleFphiProduced(PMEUnit u, bool vir)
{
   fphi_key.valid = true;
   fphi_key.vir = vir;
   fphi_key.cmp = cmp_serial;
   fphi_key.epoch = xyz_epoch;
   fphi_key.params = u->getParams();
}

bool mpoleFphiCurrent(PMEUnit u, bool need_vir)
{
   return fphi_key.valid and cmp_key.valid and fphi_key.cmp == cmp_serial and fphi_key.epoch == xyz_epoch
      and *u == fphi_key.params and (fphi_key.vir or not need_vir);
}

bool mpoleRecipVirial()
{
   return recip_vir;
}

TINKER_FVOID2(acc1, cu1, chkpole, real (*)[MPL_TOTAL], real (*)[MPL_TOTAL], real (*)[MPL_TOTAL]);
void chkpole()
{
   if (chk_epoch == xyz_epoch)
      return;
   // zaxis records one inversion per site for every copy of the multipoles, so
   // all the copies in use are inverted together, as chkpole.f does.
   auto* pole_used = (use(Potent::MPOLE) or use(Potent::POLAR)) ? pole : nullptr;
   auto* repole_used = use(Potent::REPULS) ? repole : nullptr;
   auto* poleorig_used = usePoleorig() ? poleorig : nullptr;
   TINKER_FCALL2(acc1, cu1, chkpole, pole_used, poleorig_used, repole_used);
   chk_epoch = xyz_epoch;
   // pole may have been inverted, and whatever was rotated from it with it;
   // poleorig was inverted along with it, so a scaled pole stays scaled.
   const bool scaled = pole_serial == pole_lmda_serial;
   const bool flux = flux_serial == pole_serial;
   ++pole_serial;
   if (scaled)
      pole_lmda_serial = pole_serial;
   if (flux)
      flux_serial = pole_serial;
   rpole_key = RpoleKey();
}

// The charge penetration and charge flux parameters, which altelec and altpolr
// scale along with the multipoles (mutate.f). Only free energy perturbation
// moves them after setup, a few times per sample, so they are scaled on the
// host from the originals the Fortran setup saved.
static void mpoleScaleChgpenCflux(double lmda)
{
   if (param_lmda_set and lmda == param_lmda and chargeMask() == param_lmda_mask)
      return;
   param_lmda = lmda;
   param_lmda_mask = chargeMask();
   param_lmda_set = true;

   const GrpScale sc = grpScale(lmda);
   bool wait = false;

   if (mplpot::use_chgpen) {
      std::vector<double> core(n), val(n);
      for (int i = 0; i < n; ++i) {
         double s = sc.s[mutant::mutg[i]];
         core[i] = s * dlmda::pcoreorig[i];
         val[i] = s * dlmda::pvalorig[i];
      }
      darray::copyin(g::q0, n, pcore, core.data());
      darray::copyin(g::q0, n, pval0, val.data());
      darray::copyin(g::q0, n, pval, val.data());
      wait = true;
   }

   if (use(Potent::CHGFLX)) {
      // the monopoles the flux is added to
      std::vector<double> mono(n);
      for (int i = 0; i < n; ++i)
         mono[i] = sc.s[mutant::mutg[i]] * dlmda::poleorig[mpole::maxpole * i];
      darray::copyin(g::q0, n, mono0, mono.data());

      // a bond or an angle is scaled only when the ligand holds all of it
      const int nb = bndstr::nbond, na = angbnd::nangle;
      std::vector<double> bf(nb), af(2 * na), abf(2 * na);
      for (int i = 0; i < nb; ++i) {
         int ga = mutant::mutg[bndstr::ibnd[2 * i] - 1];
         int gb = mutant::mutg[bndstr::ibnd[2 * i + 1] - 1];
         double s = (ga != 0 and gb != 0) ? sc.s[ga] : 1.0;
         bf[i] = s * dlmda::bflxorig[i];
      }
      for (int i = 0; i < na; ++i) {
         int ga = mutant::mutg[angbnd::iang[4 * i] - 1];
         int gb = mutant::mutg[angbnd::iang[4 * i + 1] - 1];
         int gc = mutant::mutg[angbnd::iang[4 * i + 2] - 1];
         double s = (ga != 0 and gb != 0 and gc != 0) ? sc.s[ga] : 1.0;
         for (int j = 0; j < 2; ++j) {
            af[2 * i + j] = s * dlmda::aflxorig[2 * i + j];
            abf[2 * i + j] = s * dlmda::abflxorig[2 * i + j];
         }
      }
      darray::copyin(g::q0, nb, bflx, bf.data());
      darray::copyin(g::q0, na, aflx, af.data());
      darray::copyin(g::q0, na, abflx, abf.data());
      wait = true;
   }

   if (wait)
      waitFor(g::q0);
}

TINKER_FVOID2(acc1, cu1, mpoleScale, const int*, GrpScale);
void mpoleScale(double lmda, bool flux_follows)
{
   // For a caller that applies the charge flux next, a pole that was scaled
   // and then given the flux is as good as the scaled one: the flux rewrites
   // the monopoles it changed. Scaling it again would only make the caller
   // redo the flux and everything built from pole since.
   bool scaled = pole_serial == pole_lmda_serial;
   if (flux_follows and pole_serial == flux_serial and flux_base_serial == pole_lmda_serial)
      scaled = true;
   // The active leg of a staged relative schedule picks the ligand that
   // charges, so the same lambda on another leg is another state.
   if (lmda == pole_lmda and chargeMask() == pole_lmda_mask and scaled)
      return;
   TINKER_FCALL2(acc1, cu1, mpoleScale, emGroup(), grpScale(lmda));
   if (use_fep)
      mpoleScaleChgpenCflux(lmda);
   pole_lmda = lmda;
   pole_lmda_mask = chargeMask();
   pole_lmda_serial = ++pole_serial;
}

static RpoleKey poleKey()
{
   RpoleKey k;
   k.src = RpoleSrc::POLE;
   k.epoch = xyz_epoch;
   k.pole = pole_serial;
   return k;
}

static RpoleKey origKey()
{
   RpoleKey k;
   k.src = RpoleSrc::ORIG;
   k.epoch = xyz_epoch;
   return k;
}

static RpoleKey stateKey(RdtMask mask, const int* group)
{
   RpoleKey k;
   k.src = RpoleSrc::MASKED;
   k.mask = mask;
   k.group = group;
   k.epoch = xyz_epoch;
   return k;
}

TINKER_FVOID2(acc1, cu1, rotpole, bool);
TINKER_FVOID2(acc1, cu1, rotpoleState, RdtMask, const int*);
static void ensureRpole(const RpoleKey& k)
{
   assert(k.src != RpoleSrc::NONE);
   if (k == rpole_key)
      return;
   if (k.src == RpoleSrc::MASKED)
      TINKER_FCALL2(acc1, cu1, rotpoleState, k.mask, k.group);
   else
      TINKER_FCALL2(acc1, cu1, rotpole, k.src == RpoleSrc::ORIG);
   rpole_key = k;
   ++rpole_serial;
}

// With do_dlmda, cmp is lambda scaled, and build_dl also fills dlcmp (d cmp / d
// lambda), which only a version carrying a lambda derivative reads.
static void ensureCmp(bool do_dlmda, bool build_dl)
{
   const CmpKey& c = cmp_key;
   const bool current = c.valid and c.rpole == rpole_serial and c.dlmda == do_dlmda
      and (not do_dlmda or (c.lmda == elam and c.group == emGroup() and (c.dl or not build_dl)));
   if (current)
      return;
   if (do_dlmda)
      rpoleToCmpDlmda(build_dl);
   else
      rpoleToCmp();
   ++cmp_serial;
   cmp_key.valid = true;
   cmp_key.dlmda = do_dlmda;
   cmp_key.dl = build_dl;
   cmp_key.lmda = do_dlmda ? elam : 0;
   cmp_key.group = do_dlmda ? emGroup() : nullptr;
   cmp_key.rpole = rpole_serial;
}

// The B-splines depend on the coordinates and the box alone, and a fill at one
// level also holds every lower level.
static void ensureSplines(PMEUnit u, int level)
{
   int iu = u;
   if (iu >= (int)spline_key.size())
      spline_key.resize(iu + 1);
   SplineKey& k = spline_key[iu];
   if (k.epoch == xyz_epoch and k.level >= level)
      return;
   bsplineFill(u, level);
   k.epoch = xyz_epoch;
   k.level = level;
}

// Only gridPut_cu2 reads precomputed B-splines: gridMpole at level 3 on every
// unit it spreads cmp on, and gridUind at level 2 on ppme_unit.
static void ensurePmeSplines(bool do_dlmda)
{
   if (not(pltfm_config & Platform::CUDA))
      return;
   constexpr bool mpole_reads = !TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE;
   constexpr bool uind_reads = !TINKER_CU_THETA_ON_THE_FLY_GRID_UIND;
   if (mpole_reads) {
      if (epme_unit.valid())
         ensureSplines(epme_unit, 3);
      // the polarization direct field spreads cmp on ppme_unit
      if (ppme_unit.valid())
         ensureSplines(ppme_unit, 3);
      if (do_dlmda and dlpme_unit.valid())
         ensureSplines(dlpme_unit, 3);
   }
   if (uind_reads and ppme_unit.valid())
      ensureSplines(ppme_unit, 2);
}

static void ensurePme(bool do_dlmda, bool build_dl)
{
   if (not useEwald())
      return;
   ensureCmp(do_dlmda, build_dl);
   ensurePmeSplines(do_dlmda);
}

void mpoleBegin(int vers, bool zero_vir_trq)
{
   recip_vir = vir_m and (vers & calc::virial);
   if (vers & calc::grad) {
      darray::zero(g::q0, n, trqx, trqy, trqz);
      if (lmdaDerivVers(vers, edlmdaActive() or pdlmdaActive()) & (calc::grad_dlmda | calc::virial_dlmda))
         darray::zero(g::q0, n, dltrqx, dltrqy, dltrqz);
   }
   if (zero_vir_trq and (vers & calc::virial))
      darray::zero(g::q0, bufferSize(), vir_trq);
   chkpole();
}

void mpoleUsePole()
{
   ensureRpole(poleKey());
   ensurePme(false, false);
}

void mpoleUseOrig(bool build_dl)
{
   ensureRpole(origKey());
   ensurePme(true, build_dl);
}

void mpoleUseState(RdtMask mask, const int* group)
{
   ensureRpole(stateKey(mask, group));
   ensurePme(false, false);
}

void mpoleInit(int vers, bool do_dlmda)
{
   mpoleBegin(vers);
   if (do_dlmda) {
      mpoleUseOrig(lmdaDerivVers(vers, edlmdaActive()) & calc::dlmda);
   } else {
      mpoleUsePole();
   }
}

void mpoleRotateScaled()
{
   ensureRpole(poleKey());
}

void mpoleEnsureElec(bool flux_follows)
{
   // A polarization state at plam leaves pole scaled there.
   if (usePoleorig())
      mpoleScale(elam, flux_follows);
}

void mpoleEnsurePhysical()
{
   mpoleEnsureElec();
   if (use(Potent::CHGFLX))
      alterchg();
   chkpole();
   ensureRpole(poleKey());
}
}
