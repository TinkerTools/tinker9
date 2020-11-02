#include "pme.h"
#include "box.h"
#include "edisp.h"
#include "elec.h"
#include "mathfunc.h"
#include "md.h"
#include "pmestuf.h"
#include "potent.h"
#include "switch.h"
#include "tinker_rt.h"
#include "tool/error.h"
#include "tool/fc.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/ewald.hh>
#include <tinker/detail/pme.hh>


namespace tinker {
bool PME::Params::operator==(const Params& st) const
{
   const double eps = 1.0e-6;
   bool ans = std::fabs(aewald - st.aewald) < eps && nfft1 == st.nfft1 &&
      nfft2 == st.nfft2 && nfft3 == st.nfft3 && bsorder == st.bsorder;
   return ans;
}


PME::Params::Params(real a, int n1, int n2, int n3, int o)
   : aewald(a)
   , nfft1(n1)
   , nfft2(n2)
   , nfft3(n3)
   , bsorder(o)
{}


void PME::set_params(const PME::Params& p)
{
   aewald = p.aewald;
   nfft1 = p.nfft1;
   nfft2 = p.nfft2;
   nfft3 = p.nfft3;
   bsorder = p.bsorder;
}


PME::Params PME::get_params() const
{
   Params p0(aewald, nfft1, nfft2, nfft3, bsorder);
   return p0;
}


bool PME::operator==(const Params& p) const
{
   return get_params() == p;
}


PME::~PME()
{
   darray::deallocate(bsmod1, bsmod2, bsmod3, qgrid);
   darray::deallocate(igrid, thetai1, thetai2, thetai3);
}

namespace {
void pme_op_alloc(PMEUnit& unit, const PME::Params& p, bool unique)
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
      darray::allocate(padded_n * p.bsorder * 4, &st.thetai1, &st.thetai2,
                       &st.thetai3);

      st.set_params(p);
   }
}


void pme_op_copyin(PMEUnit unit)
{
   if (!unit.valid())
      return;

   auto& st = *unit;

   // This code assumes that the FFT grids of an energy term will not change in
   // a calculation.
   int maxfft = max_of(st.nfft1, st.nfft2, st.nfft3);
   std::vector<double> array(st.bsorder);
   std::vector<double> bsarray(maxfft);
   double x = 0;
   TINKER_RT(bspline)(&x, &st.bsorder, array.data());
   for (int i = 0; i < maxfft; ++i) {
      bsarray[i] = 0;
   }
   assert(st.bsorder + 1 <= maxfft);
   for (int i = 0; i < st.bsorder; ++i) {
      bsarray[i + 1] = array[i];
   }
   std::vector<double> bsmodbuf(maxfft);
   TINKER_RT(dftmod)
   (bsmodbuf.data(), bsarray.data(), &st.nfft1, &st.bsorder);
   darray::copyin(async_queue, st.nfft1, st.bsmod1, bsmodbuf.data());
   wait_for(async_queue);
   TINKER_RT(dftmod)
   (bsmodbuf.data(), bsarray.data(), &st.nfft2, &st.bsorder);
   darray::copyin(async_queue, st.nfft2, st.bsmod2, bsmodbuf.data());
   wait_for(async_queue);
   TINKER_RT(dftmod)
   (bsmodbuf.data(), bsarray.data(), &st.nfft3, &st.bsorder);
   darray::copyin(async_queue, st.nfft3, st.bsmod3, bsmodbuf.data());
   wait_for(async_queue);

   unit.update_deviceptr(st, async_queue);
   wait_for(async_queue);
}
}
}


namespace tinker {
void pme_data(rc_op op)
{
   if (!use_ewald() && !use_dewald())
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
         double ecut = switch_off(switch_ewald);
         double dcut = switch_off(switch_dewald);
         double ext;
         t_extent(ext);
         double wbox = 2 * (ext + std::fmax(ecut, dcut));
         box_extent(wbox);
      }
   }


   if (use_potent(charge_term) && use_ewald()) {
      if (op & rc_alloc) {
         epme_unit.close();
         PME::Params p(ewald::aeewald, pme::nefft1, pme::nefft2, pme::nefft3,
                       pme::bseorder);
         pme_op_alloc(epme_unit, p, false);
      }


      if (op & rc_init) {
         pme_op_copyin(epme_unit);
      }
   }


   if ((use_potent(mpole_term) || use_potent(polar_term)) && use_ewald()) {
      if (op & rc_dealloc) {
         darray::deallocate(cmp, fmp, cphi, fphi);
         if (use_potent(polar_term)) {
            darray::deallocate(fuind, fuinp, fdip_phi1, fdip_phi2, cphidp,
                               fphidp);
            darray::deallocate(vir_m);
         }
      }


      if (op & rc_alloc) {
         darray::allocate(n, &cmp, &fmp, &cphi, &fphi);
         if (use_potent(polar_term)) {
            darray::allocate(n, &fuind, &fuinp, &fdip_phi1, &fdip_phi2, &cphidp,
                             &fphidp);
            if (rc_flag & calc::virial)
               darray::allocate(buffer_size(), &vir_m);
            else
               vir_m = nullptr;
         } else {
            vir_m = nullptr;
         }


         bool unique_grids = false;


         // electrostatics
         epme_unit.close();
         if (use_potent(mpole_term)) {
            unique_grids = false;
            PME::Params p(ewald::aeewald, pme::nefft1, pme::nefft2, pme::nefft3,
                          pme::bseorder);
            pme_op_alloc(epme_unit, p, unique_grids);
         }


         // polarization
         ppme_unit.close();
         pvpme_unit.close();
         if (use_potent(polar_term)) {
            PME::Params p(ewald::apewald, pme::nefft1, pme::nefft2, pme::nefft3,
                          pme::bsporder);
            pme_op_alloc(ppme_unit, p, unique_grids);
            if (rc_flag & calc::virial) {
               unique_grids = true;
               pme_op_alloc(pvpme_unit, p, unique_grids);
            }
         }
      }


      if (op & rc_init) {
         pme_op_copyin(epme_unit);
         pme_op_copyin(ppme_unit);
         pme_op_copyin(pvpme_unit);
      }
   }


   if (use_potent(disp_term) && use_dewald()) {
      if (op & rc_alloc) {
         dpme_unit.close();
         PME::Params p(ewald::adewald, pme::ndfft1, pme::ndfft2, pme::ndfft3,
                       pme::bsdorder);
         pme_op_alloc(dpme_unit, p, false);
      }


      if (op & rc_init) {
         pme_op_copyin(dpme_unit);
      }
   }


#if TINKER_CUDART
   if (op & rc_init) {
      pme_cuda_func_config();
   }
#endif
}
}
