#include "ff/energy.h"
#include "ff/hippo/cflux.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "tool/cudalib.h"
#include "tool/error.h"
#include "tool/zero.h"

#include "ff/amoeba/emplar.h"
#include "ff/amoeba/empole.h"
#include "ff/amoeba/epolar.h"
#include "ff/elec.h"
#include "ff/pchg/echarge.h"
#include "ff/pchg/echglj.h"
#include "ff/pchg/evalence.h"
#include "ff/pchg/evdw.h"
#include "ff/pme.h"

#include "ff/hippo/echgtrn.h"
#include "ff/hippo/edisp.h"
#include "ff/hippo/empolechgpen.h"
#include "ff/hippo/epolarchgpen.h"
#include "ff/hippo/erepel.h"

#include <tinker/detail/mplpot.hh>

namespace tinker {
struct DHFlow
{
   EnergyBufferTraits::type e_val;
   EnergyBufferTraits::type e_vdw;
   EnergyBufferTraits::type e_ele;
   VirialBufferTraits::type v_val[VirialBufferTraits::N];
   VirialBufferTraits::type v_vdw[VirialBufferTraits::N];
   VirialBufferTraits::type v_ele[VirialBufferTraits::N];
};

static DHFlow ev_hobj;
static DHFlow* ev_dptr;
}

namespace tinker {
static bool ecore_val;
static bool ecore_vdw;
static bool ecore_ele;

void energyData(RcOp op)
{
   if ((rc_flag & calc::vmask) == 0)
      return;

   RcMan egv42{egvData, op};

   // bonded terms

   RcMan ebond42{ebondData, op};
   RcMan eangle42{eangleData, op};
   RcMan estrbnd42{estrbndData, op};
   RcMan eurey42{eureyData, op};
   RcMan eopbend42{eopbendData, op};
   RcMan eimprop42{eimpropData, op};
   RcMan eimptor42{eimptorData, op};
   RcMan etors42{etorsData, op};
   RcMan epitors42{epitorsData, op};
   RcMan estrtor42{estrtorData, op};
   RcMan eangtor42{eangtorData, op};
   RcMan etortor42{etortorData, op};

   // misc. terms

   RcMan egeom42{egeomData, op};

   // non-bonded terms

   RcMan evdw42{evdwData, op};

   // Must call elecData() before any electrostatics routine.

   RcMan elec42{elecData, op};
   RcMan pme42{pmeData, op};

   RcMan echarge42{echargeData, op};
   // Must follow evdw_data() and echarge_data().
   RcMan echglj42{echgljData, op};

   // empoleData() must be in front of epolarData().
   RcMan empole42{empoleData, op};
   RcMan epolar42{epolarData, op};

   // HIPPO
   RcMan cflux43{cfluxData, op};
   RcMan empole43{empoleChgpenData, op};
   RcMan epolar43{epolarChgpenData, op};
   RcMan echgtrn42{echgtrnData, op};
   RcMan erepel42{erepelData, op};
   RcMan edisp42{edispData, op};

   // Must call fftData() after all of the electrostatics routines.
   RcMan fft42{fftData, op};
}

bool useEnergyVdw()
{
   bool ans = false;

   // AMOEBA
   ans = ans || usePotent(Potent::VDW);

   // HIPPO
   ans = ans || usePotent(Potent::DISP);
   ans = ans || usePotent(Potent::REPULS);

   return ans;
}

bool useEnergyElec()
{
   bool ans = false;

   // AMOEBA
   ans = ans || usePotent(Potent::CHARGE);
   ans = ans || usePotent(Potent::MPOLE);
   ans = ans || usePotent(Potent::POLAR);

   // HIPPO
   ans = ans || usePotent(Potent::CHGTRN);

   return ans;
}

const TimeScaleConfig& defaultTSConfig()
{
   static TimeScaleConfig tsconfig{
      {"ebond", 0},
      {"eangle", 0},
      {"estrbnd", 0},
      {"eurey", 0},
      {"eopbend", 0},
      {"eimprop", 0},
      {"eimptor", 0},
      {"etors", 0},
      {"epitors", 0},
      {"estrtor", 0},
      {"eangtor", 0},
      {"etortor", 0},
      {"egeom", 0},

      {"evalence", 0},

      {"evdw", 0},

      {"echarge", 0},
      {"echglj", 0},

      {"emplar", 0},
      {"empole", 0},
      {"epolar", 0},

      {"empole_chgpen", 0},
      {"epolar_chgpen", 0},

      {"echgtrn", 0},
      {"edisp", 0},
      {"erepel", 0},
      {"ehippo", 0},
   };
   return tsconfig;
}

namespace {
auto tscfg__ = [](std::string eng, bool& use_flag, unsigned tsflag,
                  const TimeScaleConfig& tsconfig) {
   auto local_flag = tsflag;
   const auto& local_cfg = tsconfig;
   try {
      bool f = local_flag & (1 << local_cfg.at(eng));
      use_flag = use_flag || f;
      return f;
   } catch (const std::out_of_range&) {
      TINKER_THROW(format("Time scale of the %s term is unknown.\n", eng));
   }
};
#define tscfg(x, f) tscfg__(x, f, tsflag, tsconfig)
}

static bool amoeba_emplar(int vers)
{
   if (mplpot::use_chgpen)
      return false;
   if (rc_flag & calc::analyz)
      return false;
   if (vers & calc::analyz)
      return false;

   return usePotent(Potent::MPOLE) && usePotent(Potent::POLAR) && (mlistVersion() & Nbl::SPATIAL);
}

static bool amoeba_empole(int vers)
{
   if (mplpot::use_chgpen)
      return false;

   if (amoeba_emplar(vers))
      return false;
   return usePotent(Potent::MPOLE);
}

static bool amoeba_epolar(int vers)
{
   if (mplpot::use_chgpen)
      return false;

   if (amoeba_emplar(vers))
      return false;
   return usePotent(Potent::POLAR);
}

static bool amoeba_echglj(int vers)
{
   if (rc_flag & calc::analyz)
      return false;
   if (vers & calc::analyz)
      return false;
   if (!usePotent(Potent::CHARGE) || !usePotent(Potent::VDW))
      return false;
   if (!(clistVersion() & Nbl::SPATIAL))
      return false;
   if (ebuffer != 0)
      return false;
   if (vdwtyp != evdw_t::lj)
      return false;
   if (vdwpr_in_use)
      return false;
   return true;
}

static bool amoeba_echarge(int vers)
{
   if (amoeba_echglj(vers))
      return false;
   return usePotent(Potent::CHARGE);
}

static bool amoeba_evdw(int vers)
{
   if (amoeba_echglj(vers))
      return false;
   return usePotent(Potent::VDW);
}

static bool hippo_empole(int vers)
{
   if (not mplpot::use_chgpen)
      return false;
   if (amoeba_emplar(vers))
      return false;
   return usePotent(Potent::MPOLE);
}

static bool hippo_epolar(int vers)
{
   if (not mplpot::use_chgpen)
      return false;
   if (amoeba_emplar(vers))
      return false;
   return usePotent(Potent::POLAR);
}

void energy_core(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   pme_stream_start_record(use_pme_stream);

   vers = vers & calc::vmask;

   ecore_val = false;
   ecore_vdw = false;
   ecore_ele = false;

   if (pltfm_config & Platform::CUDA) {
      bool calc_val = usePotent(Potent::BOND) or usePotent(Potent::ANGLE) or
         usePotent(Potent::STRBND) or usePotent(Potent::UREY) or usePotent(Potent::OPBEND) or
         usePotent(Potent::IMPROP) or usePotent(Potent::IMPTORS) or usePotent(Potent::TORSION) or
         usePotent(Potent::PITORS) or usePotent(Potent::STRTOR) or usePotent(Potent::ANGTOR) or
         usePotent(Potent::TORTOR) or usePotent(Potent::GEOM);
      if (calc_val and tscfg("evalence", ecore_val))
         evalence(vers);
   } else {
      // bonded terms

      if (usePotent(Potent::BOND))
         if (tscfg("ebond", ecore_val))
            ebond(vers);
      if (usePotent(Potent::ANGLE))
         if (tscfg("eangle", ecore_val))
            eangle(vers);
      if (usePotent(Potent::STRBND))
         if (tscfg("estrbnd", ecore_val))
            estrbnd(vers);
      if (usePotent(Potent::UREY))
         if (tscfg("eurey", ecore_val))
            eurey(vers);
      if (usePotent(Potent::OPBEND))
         if (tscfg("eopbend", ecore_val))
            eopbend(vers);
      if (usePotent(Potent::IMPROP))
         if (tscfg("eimprop", ecore_val))
            eimprop(vers);
      if (usePotent(Potent::IMPTORS))
         if (tscfg("eimptor", ecore_val))
            eimptor(vers);
      if (usePotent(Potent::TORSION))
         if (tscfg("etors", ecore_val))
            etors(vers);
      if (usePotent(Potent::PITORS))
         if (tscfg("epitors", ecore_val))
            epitors(vers);
      if (usePotent(Potent::STRTOR))
         if (tscfg("estrtor", ecore_val))
            estrtor(vers);
      if (usePotent(Potent::ANGTOR))
         if (tscfg("eangtor", ecore_val))
            eangtor(vers);
      if (usePotent(Potent::TORTOR))
         if (tscfg("etortor", ecore_val))
            etortor(vers);

      // misc. terms

      if (usePotent(Potent::GEOM))
         if (tscfg("egeom", ecore_val))
            egeom(vers);
   }

   // non-bonded terms

   if (amoeba_evdw(vers))
      if (tscfg("evdw", ecore_vdw))
         evdw(vers);

   if (amoeba_echarge(vers))
      if (tscfg("echarge", ecore_ele))
         echarge(vers);
   if (amoeba_echglj(vers))
      if (tscfg("echglj", ecore_ele)) {
         ecore_vdw = true;
         echglj(vers);
      }

   if (amoeba_empole(vers))
      if (tscfg("empole", ecore_ele))
         empole(vers);
   if (amoeba_epolar(vers))
      if (tscfg("epolar", ecore_ele))
         epolar(vers);
   if (amoeba_emplar(vers))
      if (tscfg("emplar", ecore_ele))
         emplar(vers);

   if (hippo_empole(vers))
      if (tscfg("empole_chgpen", ecore_ele))
         empoleChgpen(vers);
   if (hippo_epolar(vers))
      if (tscfg("epolar_chgpen", ecore_ele))
         epolarChgpen(vers);
   if (usePotent(Potent::CHGTRN))
      if (tscfg("echgtrn", ecore_ele))
         echgtrn(vers);
   if (usePotent(Potent::DISP))
      if (tscfg("edisp", ecore_vdw))
         edisp(vers);
   if (usePotent(Potent::REPULS))
      if (tscfg("erepel", ecore_vdw))
         erepel(vers);

   pme_stream_finish_wait(use_pme_stream and not(vers & calc::analyz));
}

void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   zero_egv(vers);
   energy_core(vers, tsflag, tsconfig);

   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   bool must_wait = false;
   ev_hobj.e_val = 0;
   ev_hobj.e_vdw = 0;
   ev_hobj.e_ele = 0;
   if (do_e) {
      if (!rc_a) {
         size_t bufsize = bufferSize();
         if (ecore_val) {
            must_wait = true;
            reduceSumOnDevice(&ev_dptr->e_val, eng_buf, bufsize, g::q0);
         }
         if (ecore_vdw && eng_buf_vdw) {
            must_wait = true;
            reduceSumOnDevice(&ev_dptr->e_vdw, eng_buf_vdw, bufsize, g::q0);
         }
         if (ecore_ele && eng_buf_elec) {
            must_wait = true;
            reduceSumOnDevice(&ev_dptr->e_ele, eng_buf_elec, bufsize, g::q0);
         }
      }
   }

   zeroOnHost(ev_hobj.v_val);
   zeroOnHost(ev_hobj.v_vdw);
   zeroOnHost(ev_hobj.v_ele);
   if (do_v) {
      if (!rc_a) {
         size_t bufsize = bufferSize();
         if (ecore_val) {
            must_wait = true;
            reduceSum2OnDevice(ev_dptr->v_val, vir_buf, bufsize, g::q0);
         }
         if (ecore_vdw && vir_buf_vdw) {
            must_wait = true;
            reduceSum2OnDevice(ev_dptr->v_vdw, vir_buf_vdw, bufsize, g::q0);
         }
         if (ecore_ele && vir_buf_elec) {
            must_wait = true;
            reduceSum2OnDevice(ev_dptr->v_ele, vir_buf_elec, bufsize, g::q0);
         }
      }
   }
   if (must_wait) {
      deviceMemoryCopyoutBytesAsync(&ev_hobj, ev_dptr, sizeof(DHFlow), g::q0);
      wait_for(g::q0);
   }
   if (do_e) {
      if (!rc_a) {
         if (ecore_val) {
            energy_valence += toFloat<energy_prec>(ev_hobj.e_val);
         }
         if (ecore_vdw && eng_buf_vdw) {
            energy_vdw += toFloat<energy_prec>(ev_hobj.e_vdw);
         }
         if (ecore_ele && eng_buf_elec) {
            energy_elec += toFloat<energy_prec>(ev_hobj.e_ele);
         }
      }
      esum = energy_valence + energy_vdw + energy_elec;
   }
   if (do_v) {
      if (!rc_a) {
         if (ecore_val) {
            virial_prec vval[VirialBufferTraits::N], v2val[9];
            for (int iv = 0; iv < (int)VirialBufferTraits::N; ++iv)
               vval[iv] = toFloat<virial_prec>(ev_hobj.v_val[iv]);
            virialReshape(v2val, vval);
            for (int iv = 0; iv < 9; ++iv)
               virial_valence[iv] += v2val[iv];
         }
         if (ecore_vdw && vir_buf_vdw) {
            virial_prec vvdw[VirialBufferTraits::N], v2vdw[9];
            for (int iv = 0; iv < (int)VirialBufferTraits::N; ++iv)
               vvdw[iv] = toFloat<virial_prec>(ev_hobj.v_vdw[iv]);
            virialReshape(v2vdw, vvdw);
            for (int iv = 0; iv < 9; ++iv)
               virial_vdw[iv] += v2vdw[iv];
         }
         if (ecore_ele && vir_buf_elec) {
            virial_prec vele[VirialBufferTraits::N], v2ele[9];
            for (int iv = 0; iv < (int)VirialBufferTraits::N; ++iv)
               vele[iv] = toFloat<virial_prec>(ev_hobj.v_ele[iv]);
            virialReshape(v2ele, vele);
            for (int iv = 0; iv < 9; ++iv)
               virial_elec[iv] += v2ele[iv];
         }
      }
      for (int iv = 0; iv < 9; ++iv)
         vir[iv] = virial_valence[iv] + virial_vdw[iv] + virial_elec[iv];
   }
   if (do_g) {
      if (ecore_vdw && gx_vdw)
         sum_gradient(gx, gy, gz, gx_vdw, gy_vdw, gz_vdw);
      if (ecore_ele && gx_elec)
         sum_gradient(gx, gy, gz, gx_elec, gy_elec, gz_elec);
   }
}

void energy(int vers)
{
   energy(vers, 1, defaultTSConfig());
}

void egvData(RcOp op)
{
   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      deviceMemoryDeallocate(ev_dptr);
      ev_dptr = nullptr;
   }

   if (op & rc_alloc) {
      deviceMemoryAllocateBytes((void**)(&ev_dptr), sizeof(DHFlow));
   }

   if (rc_flag & calc::energy) {
      if (op & rc_dealloc) {
         if (!rc_a) {
            darray::deallocate(eng_buf);
            if (useEnergyVdw())
               darray::deallocate(eng_buf_vdw);
            if (useEnergyElec())
               darray::deallocate(eng_buf_elec);
         }
      }

      if (op & rc_alloc) {
         zeroOnHost(eng_buf, eng_buf_vdw, eng_buf_elec);
         if (!rc_a) {
            auto sz = bufferSize();
            darray::allocate(sz, &eng_buf);
            if (useEnergyVdw())
               darray::allocate(sz, &eng_buf_vdw);
            if (useEnergyElec())
               darray::allocate(sz, &eng_buf_elec);
         }
      }
   }

   if (rc_flag & calc::virial) {
      if (op & rc_dealloc) {
         if (!rc_a) {
            darray::deallocate(vir_buf);
            if (useEnergyVdw())
               darray::deallocate(vir_buf_vdw);
            if (useEnergyElec())
               darray::deallocate(vir_buf_elec);
         }
      }

      if (op & rc_alloc) {
         zeroOnHost(vir_buf, vir_buf_vdw, vir_buf_elec);
         if (!rc_a) {
            auto sz = bufferSize();
            darray::allocate(sz, &vir_buf);
            if (useEnergyVdw())
               darray::allocate(sz, &vir_buf_vdw);
            if (useEnergyElec())
               darray::allocate(sz, &vir_buf_elec);
         }
      }
   }

   if (rc_flag & calc::grad) {
      if (op & rc_dealloc) {
         darray::deallocate(gx, gy, gz);
         if (useEnergyVdw())
            darray::deallocate(gx_vdw, gy_vdw, gz_vdw);
         if (useEnergyElec())
            darray::deallocate(gx_elec, gy_elec, gz_elec);
      }

      if (op & rc_alloc) {
         zeroOnHost(gx, gy, gz, gx_vdw, gy_vdw, gz_vdw, gx_elec, gy_elec, gz_elec);
         darray::allocate(n, &gx, &gy, &gz);
         if (useEnergyVdw())
            darray::allocate(n, &gx_vdw, &gy_vdw, &gz_vdw);
         if (useEnergyElec())
            darray::allocate(n, &gx_elec, &gy_elec, &gz_elec);
      }
   }
}
}
