#include "pme.h"
#include "elec.h"
#include "mathfunc.h"
#include "md.h"
#include "potent.h"
#include "tinker_rt.h"
#include <tinker/detail/ewald.hh>
#include <tinker/detail/pme.hh>

TINKER_NAMESPACE_BEGIN
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
   device_array::deallocate(bsmod1, bsmod2, bsmod3, qgrid);
   device_array::deallocate(igrid, thetai1, thetai2, thetai3);
}

static void pme_op_alloc_(PMEUnit& unit, const PME::Params& p, bool unique)
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
      device_array::allocate(p.nfft1, &st.bsmod1);
      device_array::allocate(p.nfft2, &st.bsmod2);
      device_array::allocate(p.nfft3, &st.bsmod3);
      device_array::allocate(2 * p.nfft1 * p.nfft2 * p.nfft3, &st.qgrid);
      device_array::allocate(3 * n, &st.igrid);
      device_array::allocate(padded_n * p.bsorder * 4, &st.thetai1, &st.thetai2,
                             &st.thetai3);

      st.set_params(p);
   }
}

static void pme_op_copyin_(PMEUnit unit)
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
   device_array::copyin(WAIT_NEW_Q, st.nfft1, st.bsmod1, bsmodbuf.data());
   TINKER_RT(dftmod)
   (bsmodbuf.data(), bsarray.data(), &st.nfft2, &st.bsorder);
   device_array::copyin(WAIT_NEW_Q, st.nfft2, st.bsmod2, bsmodbuf.data());
   TINKER_RT(dftmod)
   (bsmodbuf.data(), bsarray.data(), &st.nfft3, &st.bsorder);
   device_array::copyin(WAIT_NEW_Q, st.nfft3, st.bsmod3, bsmodbuf.data());

   unit.update_deviceptr(st, WAIT_NEW_Q);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void pme_init(int vers)
{
   rpole_to_cmp();

   if (vir_m)
      device_array::zero(PROCEED_NEW_Q, buffer_size(), vir_m);
}

static void pme_data1_(rc_op op)
{
   if (op & rc_dealloc) {
      PMEUnit::clear();

      device_array::deallocate(cmp, fmp, cphi, fphi);

      if (use_potent(polar_term)) {
         device_array::deallocate(fuind, fuinp, fdip_phi1, fdip_phi2, cphidp,
                                  fphidp);

         device_array::deallocate(vir_m);
      }

      epme_unit.close();
      ppme_unit.close();
      pvpme_unit.close();
      dpme_unit.close();
   }

   if (op & rc_alloc) {
      assert(PMEUnit::size() == 0);

      device_array::allocate(n, &cmp, &fmp, &cphi, &fphi);

      if (use_potent(polar_term)) {
         device_array::allocate(n, &fuind, &fuinp, &fdip_phi1, &fdip_phi2,
                                &cphidp, &fphidp);

         if (rc_flag & calc::virial) {
            device_array::allocate(buffer_size(), &vir_m);
         } else {
            vir_m = nullptr;
         }
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
         pme_op_alloc_(epme_unit, p, unique_grids);
      }

      // polarization
      ppme_unit.close();
      pvpme_unit.close();
      if (use_potent(polar_term)) {
         PME::Params p(ewald::apewald, pme::nefft1, pme::nefft2, pme::nefft3,
                       pme::bsporder);
         pme_op_alloc_(ppme_unit, p, unique_grids);
         if (rc_flag & calc::virial) {
            unique_grids = true;
            pme_op_alloc_(pvpme_unit, p, unique_grids);
         }
      }

      // dispersion
      dpme_unit.close();
      if (false) {
         unique_grids = false;
         PME::Params p(ewald::adewald, pme::ndfft1, pme::ndfft2, pme::ndfft3,
                       pme::bsdorder);
         pme_op_alloc_(dpme_unit, p, unique_grids);
      }
   }

   if (op & rc_init) {
      pme_op_copyin_(epme_unit);
      pme_op_copyin_(ppme_unit);
      pme_op_copyin_(pvpme_unit);
      pme_op_copyin_(dpme_unit);

#if TINKER_CUDART
      extern void pme_cuda_func_config();
      pme_cuda_func_config();
#endif
   }
}

void pme_data(rc_op op)
{
   if (!use_ewald())
      return;

   rc_man pme42_{pme_data1_, op};
   rc_man fft42_{fft_data, op};
}
TINKER_NAMESPACE_END
