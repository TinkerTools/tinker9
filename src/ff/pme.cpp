#include "ff/pme.h"
#include "ff/atom.h"
#include "ff/box.h"
#include "ff/elec.h"
#include "ff/hippo/edisp.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "ff/switch.h"
#include "math/maxmin.h"
#include "tool/error.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/ewald.hh>
#include <tinker/detail/pme.hh>
#include <tinker/routines.h>

namespace tinker {
bool PME::Params::operator==(const Params& st) const
{
   const double eps = 1.0e-6;
   bool ans = std::fabs(aewald - st.aewald) < eps && nfft1 == st.nfft1 && nfft2 == st.nfft2 &&
      nfft3 == st.nfft3 && bsorder == st.bsorder;
   return ans;
}

PME::Params::Params(real a, int n1, int n2, int n3, int o)
   : aewald(a)
   , nfft1(n1)
   , nfft2(n2)
   , nfft3(n3)
   , bsorder(o)
{}

void PME::setParams(const PME::Params& p)
{
   aewald = p.aewald;
   nfft1 = p.nfft1;
   nfft2 = p.nfft2;
   nfft3 = p.nfft3;
   bsorder = p.bsorder;
}

PME::Params PME::getParams() const
{
   Params p0(aewald, nfft1, nfft2, nfft3, bsorder);
   return p0;
}

bool PME::operator==(const Params& p) const
{
   return getParams() == p;
}

PME::~PME()
{
   darray::deallocate(bsmod1, bsmod2, bsmod3, qgrid);
   darray::deallocate(igrid, thetai1, thetai2, thetai3);
}

static void pmeOpAlloc(PMEUnit& unit, const PME::Params& p, bool unique)
{
   unit.close();
   for (PMEUnit idx = 0; idx < PMEUnit::size(); idx = idx + 1) {
      if (*idx == p)
         unit = idx;
   }

   if (!unit.valid() || unique == true) {
      unit = PMEUnit::open();
      auto& st = *unit;

      // see also subroutine moduli in pmestuf.f
      darray::allocate(p.nfft1, &st.bsmod1);
      darray::allocate(p.nfft2, &st.bsmod2);
      darray::allocate(p.nfft3, &st.bsmod3);
      darray::allocate(2 * p.nfft1 * p.nfft2 * p.nfft3, &st.qgrid);
      darray::allocate(3 * n, &st.igrid);
      darray::allocate(padded_n * p.bsorder * 4, &st.thetai1, &st.thetai2, &st.thetai3);

      st.setParams(p);
   }
}

static void pmeOpCopyin(PMEUnit unit)
{
   if (!unit.valid())
      return;

   auto& st = *unit;

   // This code assumes that the FFT grids of an energy term will not change in
   // a calculation.
   int maxfft = maxOf(st.nfft1, st.nfft2, st.nfft3);
   std::vector<double> array(st.bsorder);
   std::vector<double> bsarray(maxfft);
   double x = 0;
   tinker_f_bspline(&x, &st.bsorder, array.data());
   for (int i = 0; i < maxfft; ++i) {
      bsarray[i] = 0;
   }
   assert(st.bsorder + 1 <= maxfft);
   for (int i = 0; i < st.bsorder; ++i) {
      bsarray[i + 1] = array[i];
   }
   std::vector<double> bsmodbuf(maxfft);
   tinker_f_dftmod(bsmodbuf.data(), bsarray.data(), &st.nfft1, &st.bsorder);
   darray::copyin(g::q0, st.nfft1, st.bsmod1, bsmodbuf.data());
   waitFor(g::q0);
   tinker_f_dftmod(bsmodbuf.data(), bsarray.data(), &st.nfft2, &st.bsorder);
   darray::copyin(g::q0, st.nfft2, st.bsmod2, bsmodbuf.data());
   waitFor(g::q0);
   tinker_f_dftmod(bsmodbuf.data(), bsarray.data(), &st.nfft3, &st.bsorder);
   darray::copyin(g::q0, st.nfft3, st.bsmod3, bsmodbuf.data());
   waitFor(g::q0);

   unit.deviceptrUpdate(st, g::q0);
   waitFor(g::q0);
}

void pmeData(RcOp op)
{
   if (!useEwald() && !useDEwald())
      return;

   if (op & rc_dealloc) {
      PMEUnit::clear();
      epme_unit.close();
      ppme_unit.close();
      pvpme_unit.close();
      dpme_unit.close();
   }

   if (op & rc_init) {
      if (!bound::use_bounds) {
         double ecut = switchOff(Switch::EWALD);
         double dcut = switchOff(Switch::DEWALD);
         double ext;
         tinker_f_extent(&ext);
         double wbox = 2 * (ext + std::fmax(ecut, dcut));
         boxExtent(wbox);
      }
   }

   if (usePotent(Potent::CHARGE) && useEwald()) {
      if (op & rc_alloc) {
         epme_unit.close();
         PME::Params p(ewald::aeewald, pme::nefft1, pme::nefft2, pme::nefft3, pme::bseorder);
         pmeOpAlloc(epme_unit, p, false);
      }

      if (op & rc_init) {
         pmeOpCopyin(epme_unit);
      }
   }

   if ((usePotent(Potent::MPOLE) || usePotent(Potent::POLAR)) && useEwald()) {
      if (op & rc_dealloc) {
         darray::deallocate(cmp, fmp, cphi, fphi);
         if (usePotent(Potent::POLAR)) {
            darray::deallocate(fuind, fuinp, fdip_phi1, fdip_phi2, cphidp, fphidp);
            darray::deallocate(vir_m);
         }
      }

      if (op & rc_alloc) {
         darray::allocate(n, &cmp, &fmp, &cphi, &fphi);
         if (usePotent(Potent::POLAR)) {
            darray::allocate(n, &fuind, &fuinp, &fdip_phi1, &fdip_phi2, &cphidp, &fphidp);
            if (rc_flag & calc::virial)
               darray::allocate(bufferSize(), &vir_m);
            else
               vir_m = nullptr;
         } else {
            vir_m = nullptr;
         }

         bool unique_grids = false;

         // electrostatics
         epme_unit.close();
         if (usePotent(Potent::MPOLE)) {
            unique_grids = false;
            PME::Params p(ewald::aeewald, pme::nefft1, pme::nefft2, pme::nefft3, pme::bseorder);
            pmeOpAlloc(epme_unit, p, unique_grids);
         }

         // polarization
         ppme_unit.close();
         pvpme_unit.close();
         if (usePotent(Potent::POLAR)) {
            PME::Params p(ewald::apewald, pme::nefft1, pme::nefft2, pme::nefft3, pme::bsporder);
            pmeOpAlloc(ppme_unit, p, unique_grids);
            if (rc_flag & calc::virial) {
               unique_grids = true;
               pmeOpAlloc(pvpme_unit, p, unique_grids);
            }
         }
      }

      if (op & rc_init) {
         pmeOpCopyin(epme_unit);
         pmeOpCopyin(ppme_unit);
         pmeOpCopyin(pvpme_unit);
      }
   }

   if (usePotent(Potent::DISP) && useDEwald()) {
      if (op & rc_alloc) {
         dpme_unit.close();
         PME::Params p(ewald::adewald, pme::ndfft1, pme::ndfft2, pme::ndfft3, pme::bsdorder);
         pmeOpAlloc(dpme_unit, p, false);
      }

      if (op & rc_init) {
         pmeOpCopyin(dpme_unit);
      }
   }
}
}

namespace tinker {
extern void bsplineFill_cu(PMEUnit, int level);

extern void gridPchg_acc(PMEUnit, real*);
extern void gridPchg_cu(PMEUnit, real*);
extern void gridDisp_acc(PMEUnit, real*);
extern void gridDisp_cu(PMEUnit, real*);
extern void gridMpole_acc(PMEUnit, real (*)[10]);
extern void gridMpole_cu(PMEUnit, real (*)[10]);
extern void gridUind_acc(PMEUnit, real (*)[3], real (*)[3]);
extern void gridUind_cu(PMEUnit, real (*)[3], real (*)[3]);

extern void pmeConv_acc(PMEUnit, EnergyBuffer, VirialBuffer);
extern void pmeConv_cu(PMEUnit, EnergyBuffer, VirialBuffer);

extern void fphiMpole_acc(PMEUnit, real (*)[20]);
extern void fphiMpole_cu(PMEUnit, real (*)[20]);
extern void fphiUind_acc(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
extern void fphiUind_cu(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
extern void fphiUind2_acc(PMEUnit, real (*)[10], real (*)[10]);
extern void fphiUind2_cu(PMEUnit, real (*)[10], real (*)[10]);

extern void rpoleToCmp_acc();
extern void cmpToFmp_acc(PMEUnit, const real (*)[10], real (*)[10]);
extern void cuindToFuind_acc(PMEUnit, const real (*)[3], const real (*)[3], //
   real (*)[3], real (*)[3]);
extern void fphiToCphi_acc(PMEUnit, const real (*)[20], real (*)[10]);
}

namespace tinker {
void bsplineFill(PMEUnit pme_u, int level)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      bsplineFill_cu(pme_u, level);
#else
   (void)pme_u;
   (void)level;
#endif
}

void gridPchg(PMEUnit pme_u, real* pchg)
{
   int bso = pme_u->bsorder;
   if (bso != 5 and bso != 4)
      TINKER_THROW(format("gridPchg(): bsorder is %d; must be 4 or 5.\n", bso));

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      gridPchg_cu(pme_u, pchg);
   else
#endif
      gridPchg_acc(pme_u, pchg);
}

void gridMpole(PMEUnit pme_u, real (*fmp)[10])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("gridMpole(): bsorder is %d; must be 5.\n", bso));

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      gridMpole_cu(pme_u, fmp);
   else
#endif
      gridMpole_acc(pme_u, fmp);
}

void gridUind(PMEUnit pme_u, real (*fuind)[3], real (*fuinp)[3])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("gridUind(): bsorder is %d; must be 5.\n", bso));

#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      gridUind_cu(pme_u, fuind, fuinp);
   else
#endif
      gridUind_acc(pme_u, fuind, fuinp);
}

void gridDisp(PMEUnit pme_u, real* csix)
{
   int bso = pme_u->bsorder;
   if (bso != 4)
      TINKER_THROW(format("gridDisp(): bsorder is %d; must be 4.\n", bso));

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      gridDisp_cu(pme_u, csix);
   else
#endif
      gridDisp_acc(pme_u, csix);
}

void pmeConv(PMEUnit pme_u)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      pmeConv_cu(pme_u, nullptr, nullptr);
   else
#endif
      pmeConv_acc(pme_u, nullptr, nullptr);
}

void pmeConv(PMEUnit pme_u, VirialBuffer gpu_vir)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      pmeConv_cu(pme_u, nullptr, gpu_vir);
   else
#endif
      pmeConv_acc(pme_u, nullptr, gpu_vir);
}

void pmeConv(PMEUnit pme_u, EnergyBuffer gpu_e)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      pmeConv_cu(pme_u, gpu_e, nullptr);
   else
#endif
      pmeConv_acc(pme_u, gpu_e, nullptr);
}

void pmeConv(PMEUnit pme_u, EnergyBuffer gpu_e, VirialBuffer gpu_vir)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      pmeConv_cu(pme_u, gpu_e, gpu_vir);
   else
#endif
      pmeConv_acc(pme_u, gpu_e, gpu_vir);
}

void fphiMpole(PMEUnit pme_u)
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("fphiMpole(): bsorder is %d; must be 5.\n", bso));

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      fphiMpole_cu(pme_u, fphi);
   else
#endif
      fphiMpole_acc(pme_u, fphi);
}

void fphiUind(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10], real (*fdip_sum_phi)[20])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("fphiUind(): bsorder is %d; must be 5.\n", bso));

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      fphiUind_cu(pme_u, fdip_phi1, fdip_phi2, fdip_sum_phi);
   else
#endif
      fphiUind_acc(pme_u, fdip_phi1, fdip_phi2, fdip_sum_phi);
}

void fphiUind2(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("fphiUind2(): bsorder is %d; must be 5.\n", bso));

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      fphiUind2_cu(pme_u, fdip_phi1, fdip_phi2);
   else
#endif
      fphiUind2_acc(pme_u, fdip_phi1, fdip_phi2);
}

void rpoleToCmp()
{
   rpoleToCmp_acc();
}

void cmpToFmp(PMEUnit pme_u, const real (*cmp)[10], real (*fmp)[10])
{
   cmpToFmp_acc(pme_u, cmp, fmp);
}

void cuindToFuind(
   PMEUnit pme_u, const real (*cind)[3], const real (*cinp)[3], real (*fuind)[3], real (*fuinp)[3])
{
   cuindToFuind_acc(pme_u, cind, cinp, fuind, fuinp);
}

void fphiToCphi(PMEUnit pme_u, const real (*fphi)[20], real (*cphi)[10])
{
   fphiToCphi_acc(pme_u, fphi, cphi);
}
}
