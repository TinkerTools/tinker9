#include "energy.h"
#include "accmanaged.h"
#include "glob.dhflow.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/cudalib.h"
#include "tool/error.h"
#include "tool/host_zero.h"


namespace tinker {
bool ecore_val;
bool ecore_vdw;
bool ecore_ele;


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

   rc_man echarge42{echarge_data, op};
   // Must follow evdw_data() and echarge_data().
   rc_man echglj42{echglj_data, op};

   // empole_data() must be in front of epolar_data().
   rc_man empole42{empole_data, op};
   rc_man epolar42{epolar_data, op};

   // HIPPO
   rc_man echgtrn42{echgtrn_data, op};
   rc_man edisp42{edisp_data, op};

   // Must call fft_data() after all of the electrostatics routines.
   rc_man fft42{fft_data, op};
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
      {"ebond", 0},    {"eangle", 0},  {"estrbnd", 0}, {"eurey", 0},
      {"eopbend", 0},  {"eimptor", 0}, {"etors", 0},   {"epitors", 0},
      {"etortor", 0},  {"egeom", 0},

      {"evalence", 0},

      {"evdw", 0},

      {"echarge", 0},  {"echglj", 0},

      {"emplar", 0},   {"empole", 0},  {"epolar", 0},

      {"echgtrn", 0},  {"edisp", 0},   {"erepel", 0},  {"ehippo", 0},
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


void energy_core(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   pme_stream_start_record(use_pme_stream);


   vers = vers & calc::vmask;


   ecore_val = false;
   ecore_vdw = false;
   ecore_ele = false;


   if (pltfm_config & CU_PLTFM) {
      bool calc_val = use_potent(bond_term) or use_potent(angle_term) or
         use_potent(strbnd_term) or use_potent(urey_term) or
         use_potent(opbend_term) or use_potent(imptors_term) or
         use_potent(torsion_term) or use_potent(pitors_term) or
         use_potent(tortor_term) or use_potent(geom_term);
      if (calc_val and tscfg("evalence", ecore_val))
         evalence(vers);
   } else {
      // bonded terms


      if (use_potent(bond_term))
         if (tscfg("ebond", ecore_val))
            ebond(vers);
      if (use_potent(angle_term))
         if (tscfg("eangle", ecore_val))
            eangle(vers);
      if (use_potent(strbnd_term))
         if (tscfg("estrbnd", ecore_val))
            estrbnd(vers);
      if (use_potent(urey_term))
         if (tscfg("eurey", ecore_val))
            eurey(vers);
      if (use_potent(opbend_term))
         if (tscfg("eopbend", ecore_val))
            eopbend(vers);
      if (use_potent(imptors_term))
         if (tscfg("eimptor", ecore_val))
            eimptor(vers);
      if (use_potent(torsion_term))
         if (tscfg("etors", ecore_val))
            etors(vers);
      if (use_potent(pitors_term))
         if (tscfg("epitors", ecore_val))
            epitors(vers);
      if (use_potent(tortor_term))
         if (tscfg("etortor", ecore_val))
            etortor(vers);


      // misc. terms


      if (use_potent(geom_term))
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


   if (use_potent(chgtrn_term))
      if (tscfg("echgtrn", ecore_ele))
         echgtrn(vers);
   if (use_potent(disp_term))
      if (tscfg("edisp", ecore_vdw))
         edisp(vers);


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
   detail::ev_hobj.e_val = 0;
   detail::ev_hobj.e_vdw = 0;
   detail::ev_hobj.e_ele = 0;
   if (do_e) {
      if (!rc_a) {
         size_t bufsize = buffer_size();
         if (ecore_val) {
            must_wait = true;
            reduce_sum_on_device(&detail::ev_dptr->e_val, eng_buf, bufsize,
                                 g::q0);
         }
         if (ecore_vdw && eng_buf_vdw) {
            must_wait = true;
            reduce_sum_on_device(&detail::ev_dptr->e_vdw, eng_buf_vdw, bufsize,
                                 g::q0);
         }
         if (ecore_ele && eng_buf_elec) {
            must_wait = true;
            reduce_sum_on_device(&detail::ev_dptr->e_ele, eng_buf_elec, bufsize,
                                 g::q0);
         }
      }
   }


   host_zero(detail::ev_hobj.v_val);
   host_zero(detail::ev_hobj.v_vdw);
   host_zero(detail::ev_hobj.v_ele);
   if (do_v) {
      using v_t = virial_buffer_traits::type;
      if (!rc_a) {
         size_t bufsize = buffer_size();
         if (ecore_val) {
            must_wait = true;
            reduce_sum2_on_device(detail::ev_dptr->v_val, vir_buf, bufsize,
                                  g::q0);
         }
         if (ecore_vdw && vir_buf_vdw) {
            must_wait = true;
            reduce_sum2_on_device(detail::ev_dptr->v_vdw, vir_buf_vdw, bufsize,
                                  g::q0);
         }
         if (ecore_ele && vir_buf_elec) {
            must_wait = true;
            reduce_sum2_on_device(detail::ev_dptr->v_ele, vir_buf_elec, bufsize,
                                  g::q0);
         }
      }
   }
   if (must_wait) {
      device_memory_copyout_bytes_async(&detail::ev_hobj, detail::ev_dptr,
                                        sizeof(DHFlow), g::q0);
      wait_for(g::q0);
   }
   if (do_e) {
      if (!rc_a) {
         if (ecore_val) {
            energy_valence += to_flt_host<energy_prec>(detail::ev_hobj.e_val);
         }
         if (ecore_vdw && eng_buf_vdw) {
            energy_vdw += to_flt_host<energy_prec>(detail::ev_hobj.e_vdw);
         }
         if (ecore_ele && eng_buf_elec) {
            energy_elec += to_flt_host<energy_prec>(detail::ev_hobj.e_ele);
         }
      }
      esum = energy_valence + energy_vdw + energy_elec;
   }
   if (do_v) {
      if (!rc_a) {
         if (ecore_val) {
            virial_prec vval[virial_buffer_traits::N], v2val[9];
            for (int iv = 0; iv < virial_buffer_traits::N; ++iv)
               vval[iv] = to_flt_host<virial_prec>(detail::ev_hobj.v_val[iv]);
            virial_reshape(v2val, vval);
            for (int iv = 0; iv < 9; ++iv)
               virial_valence[iv] += v2val[iv];
         }
         if (ecore_vdw && vir_buf_vdw) {
            virial_prec vvdw[virial_buffer_traits::N], v2vdw[9];
            for (int iv = 0; iv < virial_buffer_traits::N; ++iv)
               vvdw[iv] = to_flt_host<virial_prec>(detail::ev_hobj.v_vdw[iv]);
            virial_reshape(v2vdw, vvdw);
            for (int iv = 0; iv < 9; ++iv)
               virial_vdw[iv] += v2vdw[iv];
         }
         if (ecore_ele && vir_buf_elec) {
            virial_prec vele[virial_buffer_traits::N], v2ele[9];
            for (int iv = 0; iv < virial_buffer_traits::N; ++iv)
               vele[iv] = to_flt_host<virial_prec>(detail::ev_hobj.v_ele[iv]);
            virial_reshape(v2ele, vele);
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
   energy(vers, 1, default_tsconfig());
}
}
