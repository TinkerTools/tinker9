#include "ff/hippo/echgtrn.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/energy.h"
#include "ff/hippo/elechippo.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/darray.h"
#include <tinker/detail/chgtrn.hh>
#include <tinker/detail/ctrpot.hh>

namespace tinker {
void echgtrnData(RcOp op)
{
   if (!usePotent(Potent::CHGTRN))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(chgct, dmpct);

      if (rc_a) {
         bufferDeallocate(calc::analyz, nct);
         bufferDeallocate(rc_flag, ect, vir_ect, dectx, decty, dectz);
      }
      nct = nullptr;
      ect = nullptr;
      vir_ect = nullptr;
      dectx = nullptr;
      decty = nullptr;
      dectz = nullptr;
   }

   if (op & rc_alloc) {
      darray::allocate(n, &chgct, &dmpct);

      nct = nullptr;
      ect = eng_buf_elec;
      vir_ect = vir_buf_elec;
      dectx = gx_elec;
      decty = gy_elec;
      dectz = gz_elec;
      if (rc_a) {
         bufferAllocate(rc_flag, &nct);
         bufferAllocate(rc_flag, &ect, &vir_ect, &dectx, &decty, &dectz);
      }
   }

   if (op & rc_init) {
      darray::copyin(g::q0, n, chgct, chgtrn::chgct);
      std::vector<real> dmpctvec(n);
      for (int i = 0; i < n; ++i) {
         real idmp = chgtrn::dmpct[i];
         if (idmp == 0)
            idmp = 1000;
         dmpctvec[i] = idmp;
      }
      darray::copyin(g::q0, n, dmpct, dmpctvec.data());
      waitFor(g::q0);
   }
}

extern void echgtrn_cu(int);
extern void echgtrn_acc(int);
void echgtrn(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   zeroOnHost(energy_ect, virial_ect);
   size_t bsize = bufferSize();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nct);
      if (do_e)
         darray::zero(g::q0, bsize, ect);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ect);
      if (do_g)
         darray::zero(g::q0, n, dectx, decty, dectz);
   }

#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      echgtrn_cu(vers);
   else
#endif
      echgtrn_acc(vers);

   if (rc_a) {
      if (do_e) {
         EnergyBuffer u = ect;
         energy_prec e = energyReduce(u);
         energy_ect += e;
         energy_elec += e;
      }
      if (do_v) {
         VirialBuffer u = vir_ect;
         virial_prec v[9];
         virialReduce(v, u);
         for (int iv = 0; iv < 9; ++iv) {
            virial_ect[iv] += v[iv];
            virial_elec[iv] += v[iv];
         }
      }
      if (do_g)
         sum_gradient(gx_elec, gy_elec, gz_elec, dectx, decty, dectz);
   }
}
}
