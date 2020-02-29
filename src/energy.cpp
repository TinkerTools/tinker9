#include "energy.h"
#include "error.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"


TINKER_NAMESPACE_BEGIN
void potential_data(rc_op op)
{
   if ((rc_flag & calc::vmask) == 0)
      return;

   rc_man egv42_{egv_data, op};

   // bonded terms

   rc_man ebond42_{ebond_data, op};
   rc_man eangle42_{eangle_data, op};
   rc_man estrbnd42_{estrbnd_data, op};
   rc_man eurey42_{eurey_data, op};
   rc_man eopbend42_{eopbend_data, op};
   rc_man etors42_{etors_data, op};
   rc_man epitors42_{epitors_data, op};
   rc_man etortor42_{etortor_data, op};

   // misc. terms

   rc_man egeom42_{egeom_data, op};

   // non-bonded terms

   rc_man evdw42_{evdw_data, op};

   // Must call elec_data() before any electrostatics routine.

   rc_man elec42_{elec_data, op};

   rc_man empole42_{empole_data, op};
   rc_man epolar42_{epolar_data, op};
#if TINKER_CUDART
   // Must follow empole_data() and epolar_data().
   rc_man emplar42_{emplar_data, op};
#endif
}


const TimeScaleConfig& default_tsconfig()
{
   static TimeScaleConfig tsconfig{
      {"ebond", 0},   {"eangle", 0},    {"estrbnd", 0},
      {"eurey", 0},   {"eopbend", 0},   {"etors", 0},
      {"epitors", 0}, {"etortor", 0},   {"egeom", 0},


      {"evdw", 0},    {"elec_init", 0}, {"torque", 0},
      {"emplar", 0},  {"empole", 0},    {"epolar", 0},
   };
   return tsconfig;
}


void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   auto TSCONFIG = [](std::string eng, int tsflag,
                      const TimeScaleConfig& tsconfig) {
      try {
         return tsflag & (1 << tsconfig.at(eng));
      } catch (const std::out_of_range&) {
         TINKER_THROW(format("Time scale of the {} term is unknown.\n", eng));
      }
   };


   vers = vers & calc::vmask;
   zero_egv(vers);


   // bonded terms


   if (use_potent(bond_term))
      if (TSCONFIG("ebond", tsflag, tsconfig))
         ebond(vers);
   if (use_potent(angle_term))
      if (TSCONFIG("eangle", tsflag, tsconfig))
         eangle(vers);
   if (use_potent(strbnd_term))
      if (TSCONFIG("estrbnd", tsflag, tsconfig))
         estrbnd(vers);
   if (use_potent(urey_term))
      if (TSCONFIG("eurey", tsflag, tsconfig))
         eurey(vers);
   if (use_potent(opbend_term))
      if (TSCONFIG("eopbend", tsflag, tsconfig))
         eopbend(vers);
   if (use_potent(torsion_term))
      if (TSCONFIG("etors", tsflag, tsconfig))
         etors(vers);
   if (use_potent(pitors_term))
      if (TSCONFIG("epitors", tsflag, tsconfig))
         epitors(vers);
   if (use_potent(tortor_term))
      if (TSCONFIG("etortor", tsflag, tsconfig))
         etortor(vers);


   // misc. terms


   if (use_potent(geom_term))
      if (TSCONFIG("egeom", tsflag, tsconfig))
         egeom(vers);


   // non-bonded terms


   if (use_potent(vdw_term))
      if (TSCONFIG("evdw", tsflag, tsconfig))
         evdw(vers);

   if (TSCONFIG("elec_init", tsflag, tsconfig))
      elec_init(vers);
#if TINKER_CUDART
   if (use_potent(mpole_term) && use_potent(polar_term) &&
       !(vers & calc::analyz) && mlist_version() == NBList::spatial) {
      if (TSCONFIG("emplar", tsflag, tsconfig))
         emplar_cu(vers);
      goto skip_mpole_polar;
   }
#endif
   if (use_potent(mpole_term))
      if (TSCONFIG("empole", tsflag, tsconfig))
         empole(vers);
   if (use_potent(polar_term))
      if (TSCONFIG("epolar", tsflag, tsconfig))
         epolar(vers);
skip_mpole_polar:
   if (TSCONFIG("torque", tsflag, tsconfig))
      torque(vers);


   sum_energies(vers);
}


void energy(int vers)
{
   energy(vers, 1, default_tsconfig());
}


void copy_energy(int vers, energy_prec* restrict eng, real* restrict grdx,
                 real* restrict grdy, real* restrict grdz,
                 real* restrict virial)
{
   if (eng && vers & calc::energy && eng != &esum) {
      eng[0] = esum;
   }


   if (vers & calc::grad) {
      if (grdx && grdx != gx)
         device_array::copy(PROCEED_NEW_Q, n, grdx, gx);
      if (grdy && grdy != gy)
         device_array::copy(PROCEED_NEW_Q, n, grdy, gy);
      if (grdz && grdz != gz)
         device_array::copy(PROCEED_NEW_Q, n, grdz, gz);
   }


   if (virial && vers & calc::virial && virial != &vir[0]) {
      for (int i = 0; i < 9; ++i)
         virial[i] = vir[i];
   }
}
TINKER_NAMESPACE_END
