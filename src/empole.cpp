#include "ff/amoeba/empole.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "math/zero.h"
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/sizes.hh>

namespace tinker {
void empole_nonewald_acc(int vers);
void empole_ewald_recip_acc(int vers);
void empole_ewald_real_self_acc(int vers);
void empole_nonewald_cu(int vers);
void empole_ewald_real_self_cu(int vers);

void empole_nonewald(int vers);
void empole_ewald(int vers);
void empole_ewald_real_self(int vers);
void empole_ewald_recip(int vers);
}

namespace tinker {
void empoleData(RcOp op)
{
   if (!usePotent(Potent::MPOLE))
      return;
   if (mplpot::use_chgpen)
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      if (rc_a) {
         bufferDeallocate(rc_flag, nem);
         bufferDeallocate(rc_flag, em, vir_em, demx, demy, demz);
      }
      nem = nullptr;
      em = nullptr;
      vir_em = nullptr;
      demx = nullptr;
      demy = nullptr;
      demz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      nem = nullptr;
      em = eng_buf_elec;
      vir_em = vir_buf_elec;
      demx = gx_elec;
      demy = gy_elec;
      demz = gz_elec;
      if (rc_a) {
         bufferAllocate(rc_flag, &nem);
         bufferAllocate(rc_flag, &em, &vir_em, &demx, &demy, &demz);
      }
   }

   if (op & RcOp::INIT) {}
}

void empole(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   zeroOnHost(energy_em, virial_em);
   size_t bsize = bufferSize();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nem);
      if (do_e)
         darray::zero(g::q0, bsize, em);
      if (do_v)
         darray::zero(g::q0, bsize, vir_em);
      if (do_g)
         darray::zero(g::q0, n, demx, demy, demz);
   }

   mpoleInit(vers);
   if (useEwald())
      empole_ewald(vers);
   else
      empole_nonewald(vers);
   torque(vers, demx, demy, demz);
   if (do_v) {
      VirialBuffer u2 = vir_trq;
      virial_prec v2[9];
      virialReduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv) {
         virial_em[iv] += v2[iv];
         virial_elec[iv] += v2[iv];
      }
   }

   if (rc_a) {
      if (do_e) {
         EnergyBuffer u = em;
         energy_prec e = energyReduce(u);
         energy_em += e;
         energy_elec += e;
      }
      if (do_v) {
         VirialBuffer u1 = vir_em;
         virial_prec v1[9];
         virialReduce(v1, u1);
         for (int iv = 0; iv < 9; ++iv) {
            virial_em[iv] += v1[iv];
            virial_elec[iv] += v1[iv];
         }
      }
      if (do_g)
         sumGradient(gx_elec, gy_elec, gz_elec, demx, demy, demz);
   }
}

void empole_nonewald(int vers)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      empole_nonewald_cu(vers);
   else
#endif
      empole_nonewald_acc(vers);
}

void empole_ewald(int vers)
{
   empole_ewald_real_self(vers);
   empole_ewald_recip(vers);
}

void empole_ewald_real_self(int vers)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      empole_ewald_real_self_cu(vers);
   else
#endif
      empole_ewald_real_self_acc(vers);
}

void empole_ewald_recip(int vers)
{
   empole_ewald_recip_acc(vers);
}
}
