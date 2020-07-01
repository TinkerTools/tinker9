#include "energy.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/error.h"


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
   rc_man eimptor42{eimptor_data, op};
   rc_man etors42{etors_data, op};
   rc_man epitors42{epitors_data, op};
   rc_man etortor42{etortor_data, op};

   // misc. terms

   rc_man egeom42{egeom_data, op};

   // non-bonded terms

   rc_man evdw42{evdw_data, op};

   // Must call elec_data() before any electrostatics routine.

   rc_man elec42{elec_data, op};
   rc_man pme42{pme_data, op};
   rc_man fft42{fft_data, op};

   rc_man echarge42{echarge_data, op};

   // empole_data() must be in front of epolar_data().
   rc_man empole42{empole_data, op};
   rc_man epolar42{epolar_data, op};
   // Must follow empole_data() and epolar_data().
   rc_man emplar42{emplar_data, op};

   // HIPPO
   rc_man echgtrn42{echgtrn_data, op};
   rc_man erepel42{erepel_data, op};
   rc_man edisp42{edisp_data, op};
}


bool use_energi_vdw()
{
   bool ans = false;

   // AMOEBA
   ans = ans || use_potent(vdw_term);

   // HIPPO
   ans = ans || use_potent(disp_term);
   ans = ans || use_potent(repuls_term);

   return ans;
}


bool use_energi_elec()
{
   bool ans = false;

   // AMOEBA
   ans = ans || use_potent(charge_term);
   ans = ans || use_potent(mpole_term);
   ans = ans || use_potent(polar_term);

   // HIPPO
   ans = ans || use_potent(chgtrn_term);

   return ans;
}


const TimeScaleConfig& default_tsconfig()
{
   static TimeScaleConfig tsconfig{
      {"ebond", 0},   {"eangle", 0},  {"estrbnd", 0}, {"eurey", 0},
      {"eopbend", 0}, {"eimptor", 0}, {"etors", 0},   {"epitors", 0},
      {"etortor", 0}, {"egeom", 0},

      {"evdw", 0},

      {"echarge", 0},

      {"emplar", 0},  {"empole", 0},  {"epolar", 0},

      {"echgtrn", 0}, {"edisp", 0},   {"erepel", 0},  {"ehippo", 0},
   };
   return tsconfig;
}


namespace {
auto tscfg__ = [](std::string eng, unsigned tsflag,
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
   if (use_potent(imptors_term))
      if (tscfg("eimptor"))
         eimptor(vers);
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


   if (amoeba_empole(vers))
      if (tscfg("empole"))
         empole(vers);
   if (amoeba_epolar(vers))
      if (tscfg("epolar"))
         epolar(vers);
   if (amoeba_emplar(vers))
      if (tscfg("emplar"))
         emplar(vers);


   if (use_potent(chgtrn_term))
      if (tscfg("echgtrn"))
         echgtrn(vers);
   if (use_potent(repuls_term))
      if (tscfg("erepel"))
         erepel(vers);
   if (use_potent(disp_term))
      if (tscfg("edisp"))
         edisp(vers);
}


void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   zero_egv(vers);
   energy_core(vers, tsflag, tsconfig);


   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (do_e) {
      if (!rc_a) {
         energy_prec e;
         e = energy_reduce(eng_buf);
         energy_valence += e;
         if (eng_buf_vdw) {
            e = energy_reduce(eng_buf_vdw);
            energy_vdw += e;
         }
         if (eng_buf_elec) {
            e = energy_reduce(eng_buf_elec);
            energy_elec += e;
         }
      }
      esum = energy_valence + energy_vdw + energy_elec;
   }


   if (do_v) {
      if (!rc_a) {
         virial_prec v[9];
         virial_reduce(v, vir_buf);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += v[iv];
         if (vir_buf_vdw) {
            virial_reduce(v, vir_buf_vdw);
            for (int iv = 0; iv < 9; ++iv)
               virial_vdw[iv] += v[iv];
         }
         if (vir_buf_elec) {
            virial_reduce(v, vir_buf_elec);
            for (int iv = 0; iv < 9; ++iv)
               virial_elec[iv] += v[iv];
         }
      }
      for (int iv = 0; iv < 9; ++iv)
         vir[iv] = virial_valence[iv] + virial_vdw[iv] + virial_elec[iv];
   }


   if (do_g) {
      if (gx_vdw)
         sum_gradient(gx, gy, gz, gx_vdw, gy_vdw, gz_vdw);
      if (gx_elec)
         sum_gradient(gx, gy, gz, gx_elec, gy_elec, gz_elec);
   }
}


void energy(int vers)
{
   energy(vers, 1, default_tsconfig());
}
}
