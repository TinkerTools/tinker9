#include "energy.h"
#include "error.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"


namespace tinker {
void energy_data(rc_op op)
{
   if ((rc_flag & calc::vmask) == 0)
      return;

   rc_man egv42{egv_data, op};

   // bonded terms

   rc_man ebond42{ebond_data, op};
   rc_man eangle42{eangle_data, op};
   rc_man estrbnd42{estrbnd_data, op};
   rc_man eurey42{eurey_data, op};
   rc_man eopbend42{eopbend_data, op};
   rc_man etors42{etors_data, op};
   rc_man epitors42{epitors_data, op};
   rc_man etortor42{etortor_data, op};

   // misc. terms

   rc_man egeom42{egeom_data, op};

   // non-bonded terms

   rc_man evdw42{evdw_data, op};

   // Must call elec_data() before any electrostatics routine.

   rc_man elec42{elec_data, op};

   rc_man echarge42{echarge_data, op};

   // empole_data() must be in front of epolar_data().
   rc_man empole42{empole_data, op};
   rc_man epolar42{epolar_data, op};
   // Must follow empole_data() and epolar_data().
   rc_man emplar42{emplar_data, op};
}


const TimeScaleConfig& default_tsconfig()
{
   static TimeScaleConfig tsconfig{
      {"ebond", 0},   {"eangle", 0},  {"estrbnd", 0},
      {"eurey", 0},   {"eopbend", 0}, {"etors", 0},
      {"epitors", 0}, {"etortor", 0}, {"egeom", 0},

      {"evdw", 0},

      {"echarge", 0},

      {"emplar", 0},  {"empole", 0},  {"epolar", 0},
   };
   return tsconfig;
}


namespace {
auto tscfg__ = [&](std::string eng, unsigned tsflag,
                   const TimeScaleConfig& tsconfig) {
   auto local_flag = tsflag;
   const auto& local_cfg = tsconfig;
   try {
      return local_flag & (1 << local_cfg.at(eng));
   } catch (const std::out_of_range&) {
      TINKER_THROW(format("Time scale of the %s term is unknown.\n", eng));
   }
};
#define tscfg(x) tscfg__(x, tsflag, tsconfig)
}


void energy_core(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   vers = vers & calc::vmask;


   // bonded terms


   if (use_potent(bond_term))
      if (tscfg("ebond"))
         ebond(vers);
   if (use_potent(angle_term))
      if (tscfg("eangle"))
         eangle(vers);
   if (use_potent(strbnd_term))
      if (tscfg("estrbnd"))
         estrbnd(vers);
   if (use_potent(urey_term))
      if (tscfg("eurey"))
         eurey(vers);
   if (use_potent(opbend_term))
      if (tscfg("eopbend"))
         eopbend(vers);
   if (use_potent(torsion_term))
      if (tscfg("etors"))
         etors(vers);
   if (use_potent(pitors_term))
      if (tscfg("epitors"))
         epitors(vers);
   if (use_potent(tortor_term))
      if (tscfg("etortor"))
         etortor(vers);


   // misc. terms


   if (use_potent(geom_term))
      if (tscfg("egeom"))
         egeom(vers);


   // non-bonded terms


   if (use_potent(vdw_term))
      if (tscfg("evdw"))
         evdw(vers);


   if (use_potent(charge_term))
      if (tscfg("echarge"))
         echarge(vers);


   bool calc_mpole = use_potent(mpole_term); // quadrupole
   bool calc_polar = use_potent(polar_term); // AMOEBA polarization (Thole)
   bool calc_mplar = calc_mpole && calc_polar;
   calc_mplar = calc_mplar && !(vers & calc::analyz);
   calc_mplar = calc_mplar && (mlist_version() & NBL_SPATIAL);
   if (calc_mplar) {
      if (tscfg("emplar"))
         emplar(vers);
   } else {
      // update logical flags by time-scale configurations
      calc_mpole = calc_mpole && tscfg("empole");
      calc_polar = calc_polar && tscfg("epolar");
      if (calc_mpole)
         empole(vers);
      if (calc_polar)
         epolar(vers);
   }
}


void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   zero_egv(vers);
   energy_core(vers, tsflag, tsconfig);


   bool rc_a = rc_flag & calc::analyz;
   if (vers & calc::energy) {
      if (!rc_a) {
         energy_valence = energy_reduce(eng_buf);
      }
      esum = energy_valence + energy_ev + energy_elec;
   }


   if (vers & calc::virial) {
      if (!rc_a) {
         virial_reduce(virial_valence, vir_buf);
      }
      for (int iv = 0; iv < 9; ++iv)
         vir[iv] = virial_valence[iv] + virial_ev[iv] + virial_elec[iv];
   }


   if (vers & calc::grad) {
      if (devx && devy && devz)
         sum_gradient(gx, gy, gz, devx, devy, devz);
      if (decx && decy && decz)
         sum_gradient(gx, gy, gz, decx, decy, decz);
      if (demx && demy && demz)
         sum_gradient(gx, gy, gz, demx, demy, demz);
      if (depx && depy && depz && (depx != demx) && (depy != demy) &&
          (depz != demz))
         sum_gradient(gx, gy, gz, depx, depy, depz);
   }
}


void energy(int vers)
{
   energy(vers, 1, default_tsconfig());
}
}
