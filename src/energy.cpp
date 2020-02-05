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
   static TimeScaleConfig tsconfig;
   static bool init = false;
   if (!init) {
      tsconfig["ebond"] = 0;
      tsconfig["eangle"] = 0;
      tsconfig["estrbnd"] = 0;
      tsconfig["eurey"] = 0;
      tsconfig["eopbend"] = 0;
      tsconfig["etors"] = 0;
      tsconfig["epitors"] = 0;
      tsconfig["etortor"] = 0;


      tsconfig["egeom"] = 0;


      tsconfig["evdw"] = 0;


      tsconfig["elec_init"] = 0;
      tsconfig["torque"] = 0;
      tsconfig["emplar"] = 0;
      tsconfig["empole"] = 0;
      tsconfig["epolar"] = 0;


      init = true;
   }
   return tsconfig;
}


void energy_potential(int vers, int time_scale, const TimeScaleConfig& tsconfig)
{
   auto TSCONFIG = [&](const char* eng) {
      try {
         return tsconfig.at(eng);
      } catch (const std::out_of_range&) {
         TINKER_THROW(format("Time scale of the {} term is unknown.\n", eng));
      }
   };


   zero_egv(vers);


   // bonded terms


   if (use_potent(bond_term))
      if (time_scale & (1 << TSCONFIG("ebond")))
         ebond(vers);
   if (use_potent(angle_term))
      if (time_scale & (1 << TSCONFIG("eangle")))
         eangle(vers);
   if (use_potent(strbnd_term))
      if (time_scale & (1 << TSCONFIG("estrbnd")))
         estrbnd(vers);
   if (use_potent(urey_term))
      if (time_scale & (1 << TSCONFIG("eurey")))
         eurey(vers);
   if (use_potent(opbend_term))
      if (time_scale & (1 << TSCONFIG("eopbend")))
         eopbend(vers);
   if (use_potent(torsion_term))
      if (time_scale & (1 << TSCONFIG("etors")))
         etors(vers);
   if (use_potent(pitors_term))
      if (time_scale & (1 << TSCONFIG("epitors")))
         epitors(vers);
   if (use_potent(tortor_term))
      if (time_scale & (1 << TSCONFIG("etortor")))
         etortor(vers);


   // misc. terms


   if (use_potent(geom_term))
      if (time_scale & (1 << TSCONFIG("egeom")))
         egeom(vers);


   // non-bonded terms


   if (use_potent(vdw_term))
      if (time_scale & (1 << TSCONFIG("evdw")))
         evdw(vers);


   elec_init(vers);
#if TINKER_CUDART
   if (use_potent(mpole_term) && use_potent(polar_term) &&
       !(vers & calc::analyz) && mlist_version() == NBList::spatial) {
      if (time_scale & (1 << TSCONFIG("emplar")))
         emplar_cu(vers);
      goto skip_mpole_polar;
   }
#endif
   if (use_potent(mpole_term))
      if (time_scale & (1 << TSCONFIG("empole")))
         empole(vers);
   if (use_potent(polar_term))
      if (time_scale & (1 << TSCONFIG("epolar")))
         epolar(vers);
skip_mpole_polar:
   torque(vers);


   sum_energies(vers);
}
TINKER_NAMESPACE_END
