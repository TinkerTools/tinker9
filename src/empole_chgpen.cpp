#include "ff/amoeba/empole.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/hippo/cflux.h"
#include "ff/hippo/empolechgpen.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "ff/amoeba/elecamoeba.h"
#include "tool/zero.h"
#include <tinker/detail/chgpen.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/potent.hh>
#include <tinker/detail/sizes.hh>

namespace tinker {
void empole_chgpen_nonewald_acc(int vers, int use_cf);
void empole_chgpen_ewald_recip_acc(int vers, int use_cf);
void empole_chgpen_ewald_real_self_acc(int vers, int use_cf);
void empole_chgpen_nonewald_cu(int vers, int use_cf);
void empole_chgpen_ewald_real_self_cu(int vers, int use_cf);

void empole_chgpen_nonewald(int vers, int use_cf);
void empole_chgpen_ewald(int vers, int use_cf);
void empole_chgpen_ewald_real_self(int vers, int use_cf);
void empole_chgpen_ewald_recip(int vers, int use_cf);
}

namespace tinker {
void empoleChgpenData(RcOp op)
{
   if (not usePotent(Potent::MPOLE) and not usePotent(Potent::CHGTRN))
      return;
   if (not mplpot::use_chgpen)
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
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

   if (op & rc_alloc) {
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

   if (op & rc_init) {}
}

void empoleChgpen(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;
   int use_cf = potent::use_chgflx;
   int use_cfgrad = use_cf and do_g;

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

   if (use_cf)
      alterchg();
   mpoleInit(vers);
   if (use_cfgrad) {
      cfluxZeroPot();
   }
   if (useEwald())
      empole_chgpen_ewald(vers, use_cfgrad);
   else
      empole_chgpen_nonewald(vers, use_cfgrad);
   torque(vers, demx, demy, demz);
   if (use_cfgrad)
      dcflux(vers, demx, demy, demz, vir_em);
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
         sum_gradient(gx_elec, gy_elec, gz_elec, demx, demy, demz);
   }
}

void empole_chgpen_nonewald(int vers, int use_cf)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      empole_chgpen_nonewald_cu(vers, use_cf);
   else
#endif
      empole_chgpen_nonewald_acc(vers, use_cf);
}

void empole_chgpen_ewald(int vers, int use_cf)
{
   empole_chgpen_ewald_real_self(vers, use_cf);
   empole_chgpen_ewald_recip(vers, use_cf);
}

void empole_chgpen_ewald_real_self(int vers, int use_cf)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      empole_chgpen_ewald_real_self_cu(vers, use_cf);
   else
#endif
      empole_chgpen_ewald_real_self_acc(vers, use_cf);
}

void empole_chgpen_ewald_recip(int vers, int use_cf)
{
   empole_chgpen_ewald_recip_acc(vers, use_cf);
}
}
