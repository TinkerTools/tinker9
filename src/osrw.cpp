#include "osrw.h"
#include "box.h"
#include "energy.h"
#include "evdw.h"
#include "fc.h"
#include "potent.h"
#include "tinker_rt.h"


namespace tinker {
bool use_osrw;
double osrw_lambda;
int osrw_vdw;
int osrw_ele;
int osrw_tor;


energy_prec osrw_du1;
virial_prec osrw_dv1[9];
grad_prec *osrw_dgx, *osrw_dgy, *osrw_dgz;


real* osrw_pchg;
real (*osrw_pole)[mpl_total];
real* osrw_polarity;
int osrw_ntbnd;
int (*osrw_itbnd)[2];
real (*osrw_tors1)[4];
real (*osrw_tors2)[4];
real (*osrw_tors3)[4];
real (*osrw_tors4)[4];
real (*osrw_tors5)[4];
real (*osrw_tors6)[4];


namespace {
grad_prec *osrw_gx, *osrw_gy, *osrw_gz;
}


void osrw_mech()
{
   use_osrw = false;
   osrw_lambda = 1.0;
   osrw_vdw = OSRW_LAM_LINEAR;
   osrw_ele = OSRW_LAM_LINEAR;
   osrw_tor = OSRW_LAM_LINEAR;


   osrw_du1 = 0;
   for (int iv = 0; iv < 9; ++iv)
      osrw_dv1[iv] = 0;
   osrw_dgx = nullptr;
   osrw_dgy = nullptr;
   osrw_dgz = nullptr;


   osrw_pchg = nullptr;
   osrw_pole = nullptr;
   osrw_polarity = nullptr;
   osrw_ntbnd = 0;
   osrw_itbnd = nullptr;
   osrw_tors1 = nullptr;
   osrw_tors2 = nullptr;
   osrw_tors3 = nullptr;
   osrw_tors4 = nullptr;
   osrw_tors5 = nullptr;
   osrw_tors6 = nullptr;


   osrw_gx = nullptr;
   osrw_gy = nullptr;
   osrw_gz = nullptr;


   // OSRW keywords.


   get_kbool("OSRW-LAMBDA", use_osrw, false);
   get_kv("OSRW-LAMBDA", osrw_lambda, 1.0);


   auto assign = [](int& func_int, const std::string& name) {
      if (name == "LINEAR")
         func_int = OSRW_LAM_LINEAR;
      else if (name == "QUADRATIC")
         func_int = OSRW_LAM_QUADRATIC;
   };
   std::string record;
   get_kv("OSRW-VDW", record, "LINEAR");
   assign(osrw_vdw, record);
   get_kv("OSRW-ELE", record, "LINEAR");
   assign(osrw_ele, record);
   get_kv("OSRW-TORS", record, "LINEAR");
   assign(osrw_tor, record);
}


void osrw_data(rc_op op)
{
   if (!use_osrw)
      return;


   if (op & rc_dealloc) {
      if (rc_flag & calc::grad) {
         darray::deallocate(osrw_dgx, osrw_dgy, osrw_dgz);
         darray::deallocate(osrw_gx, osrw_gy, osrw_gz);
      }


      if (use_potent(torsion_term)) {
         if (osrw_ntbnd > 0) {
            darray::deallocate(osrw_itbnd, osrw_tors1, osrw_tors2, osrw_tors3,
                               osrw_tors4, osrw_tors5, osrw_tors6);
            osrw_ntbnd = 0;
         }
      }


      if (use_potent(charge_term)) {
         darray::deallocate(osrw_pchg);
      }


      if (use_potent(mpole_term) || use_potent(polar_term)) {
         darray::deallocate(osrw_pole);
         if (use_potent(polar_term)) {
            darray::deallocate(osrw_polarity);
         }
      }
   }


   if (op & rc_alloc) {
      if (rc_flag & calc::grad) {
         darray::allocate(n, &osrw_dgx, &osrw_dgy, &osrw_dgz);
         darray::allocate(n, &osrw_gx, &osrw_gy, &osrw_gz);
      }


      if (use_potent(torsion_term)) {
         std::vector<int> buf;
         std::vector<std::string> vs;
         get_kv("ROTATABLE-BOND", vs, "");
         for (size_t i = 0; i < vs.size(); i += 2) {
            std::string s1 = vs.at(i);
            std::string s2 = vs.at(i + 1);
            buf.push_back(std::stoi(s1) - 1);
            buf.push_back(std::stoi(s2) - 1);
            osrw_ntbnd += 1;
         }
         if (osrw_ntbnd > 0) {
            darray::allocate(osrw_ntbnd, &osrw_itbnd);
            darray::allocate(ntors, &osrw_tors1, &osrw_tors2, &osrw_tors3,
                             &osrw_tors4, &osrw_tors5, &osrw_tors6);
            darray::copy(PROCEED_NEW_Q, ntors, osrw_tors1, tors1);
            darray::copy(PROCEED_NEW_Q, ntors, osrw_tors2, tors2);
            darray::copy(PROCEED_NEW_Q, ntors, osrw_tors3, tors3);
            darray::copy(PROCEED_NEW_Q, ntors, osrw_tors4, tors4);
            darray::copy(PROCEED_NEW_Q, ntors, osrw_tors5, tors5);
            darray::copy(PROCEED_NEW_Q, ntors, osrw_tors6, tors6);
            darray::copyin(WAIT_NEW_Q, osrw_ntbnd, osrw_itbnd, buf.data());
         }
      }


      if (use_potent(charge_term)) {
         darray::allocate(n, &osrw_pchg);
      }


      if (use_potent(mpole_term) || use_potent(polar_term)) {
         darray::allocate(n, &osrw_pole);
         if (use_potent(polar_term)) {
            darray::allocate(n, &osrw_polarity);
         }
      }
   }


   if (op & rc_init) {
      if (use_potent(charge_term)) {
         darray::copy(PROCEED_NEW_Q, n, osrw_pchg, pchg);
      }


      if (use_potent(mpole_term) || use_potent(polar_term)) {
         darray::copy(PROCEED_NEW_Q, n, osrw_pole, pole);
         if (use_potent(polar_term)) {
            darray::copy(PROCEED_NEW_Q, n, osrw_polarity, polarity);
         }
      }
   }
}


double osrw_lam_expr0(int form, double lam)
{
   if (lam <= 0)
      return 0;
   else if (lam >= 1)
      return 1;


   double ans;
   switch (form) {
   case OSRW_LAM_QUADRATIC:
      ans = lam * lam;
      break;
   default: // LINEAR
      ans = lam;
      break;
   }
   return ans;
}


double osrw_lam_expr1(int form, double lam)
{
   double ans;
   switch (form) {
   case OSRW_LAM_QUADRATIC:
      ans = 2 * lam;
      break;
   default: // LINEAR
      ans = 1;
      break;
   }
   return ans;
}


double osrw_lam_expr2(int form, double lam)
{
   double ans;
   switch (form) {
   case OSRW_LAM_QUADRATIC:
      ans = 2;
      break;
   default: // LINEAR
      ans = 0;
      break;
   }
   return ans;
}


void osrw_altele(double el)
{
   osrw_altele_acc(el);
}


void osrw_alttor(double tl)
{
   osrw_alttor_acc(tl);
}


void osrw_altvdw(double vl)
{
   vlam = vl;
}


void osrw_energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig)
{
   double sele1 = osrw_lam_expr1(osrw_ele, osrw_lambda);
   double svdw1 = osrw_lam_expr1(osrw_vdw, osrw_lambda);
   double stor1 = osrw_lam_expr1(osrw_tor, osrw_lambda);


   energy_prec osrw_du0 = 0;
   osrw_du1 = 0;
   virial_prec osrw_dv0[9];
   for (int iv = 0; iv < 9; ++iv) {
      osrw_dv0[iv] = 0;
      osrw_dv1[iv] = 0;
   }
   if (vers & calc::grad) {
      darray::zero(PROCEED_NEW_Q, n, osrw_gx, osrw_gy, osrw_gz);
      darray::zero(PROCEED_NEW_Q, n, osrw_dgx, osrw_dgy, osrw_dgz);
   }


   double sele0;
   double svdw0;
   double stor0;


   sele0 = 1 - osrw_lam_expr0(osrw_ele, osrw_lambda);
   svdw0 = 1 - osrw_lam_expr0(osrw_vdw, osrw_lambda);
   stor0 = 1 - osrw_lam_expr0(osrw_tor, osrw_lambda);
   osrw_alttor(0);
   osrw_altvdw(0);
   osrw_altele(0);


   zero_egv(vers);
   energy_core(vers, tsflag, tsconfig);
   if (vers & calc::energy) {
      for (size_t i = 0; i < energy_buffers.size(); ++i) {
         energy_buffer u = energy_buffers[i];
         energy_prec e = energy_reduce(u);
         energy_prec* eptr = get_energy_reduce_dst(u);
         double scal0 = stor0;
         double scal1 = stor1;
         if (eptr == &energy_ev) {
            scal0 = svdw0;
            scal1 = svdw1;
         } else if (eptr == &energy_ec || eptr == &energy_em ||
                    eptr == &energy_ep) {
            scal0 = sele0;
            scal1 = sele1;
         }
         osrw_du0 += scal0 * e;
         osrw_du1 -= scal1 * e;


         e *= scal0;
         *eptr = e;
      }
   }
   if (vers & calc::virial) {
      for (size_t i = 0; i < virial_buffers.size(); ++i) {
         virial_buffer u = virial_buffers[i];
         virial_prec v[9];
         virial_reduce(v, u);
         double scal0 = stor0;
         double scal1 = stor1;
         if (u == vir_ev) {
            scal0 = svdw0;
            scal1 = svdw1;
         } else if (u == vir_ec || u == vir_em || u == vir_ep) {
            scal0 = sele0;
            scal1 = sele1;
         }
         for (int iv = 0; iv < 9; ++iv) {
            osrw_dv0[iv] += scal0 * v[iv];
            osrw_dv1[iv] -= scal1 * v[iv];
         }
      }
      for (int iv = 0; iv < 9; ++iv) {
         vir[iv] = osrw_dv0[iv];
      }
   }
   if (vers & calc::grad) {
      double scale = stor0;
      scale_gradient(scale, gx, gy, gz);
      size_t ngrad = x_grads.size();
      for (size_t i = 1; i < ngrad; ++i) {
         grad_prec* gx1 = x_grads[i];
         grad_prec* gy1 = y_grads[i];
         grad_prec* gz1 = z_grads[i];
         double scal0 = stor0;
         double scal1 = stor1;
         if (gx1 == devx) {
            scal0 = svdw0;
            scal1 = svdw1;
         } else if (gx1 == decx || gx1 == demx || gx1 == depx) {
            scal0 = sele0;
            scal1 = sele1;
         }
         sum_gradient(-scal1, osrw_dgx, osrw_dgy, osrw_dgz, gx1, gy1, gz1);
         scale_gradient(scal0, gx1, gy1, gz1);
         sum_gradient(gx, gy, gz, gx1, gy1, gz1);
      }
      darray::copy(PROCEED_NEW_Q, n, osrw_gx, gx);
      darray::copy(PROCEED_NEW_Q, n, osrw_gy, gy);
      darray::copy(PROCEED_NEW_Q, n, osrw_gz, gz);
   }


   sele0 = 1.0 - sele0;
   svdw0 = 1.0 - svdw0;
   stor0 = 1.0 - stor0;
   osrw_alttor(1);
   osrw_altvdw(1);
   osrw_altele(1);


   zero_egv(vers);
   energy_core(vers, tsflag, tsconfig);
   if (vers & calc::energy) {
      for (size_t i = 0; i < energy_buffers.size(); ++i) {
         energy_buffer u = energy_buffers[i];
         energy_prec e = energy_reduce(u);
         energy_prec* eptr = get_energy_reduce_dst(u);
         double scal0 = stor0;
         double scal1 = stor1;
         if (eptr == &energy_ev) {
            scal0 = svdw0;
            scal1 = svdw1;
         } else if (eptr == &energy_ec || eptr == &energy_em ||
                    eptr == &energy_ep) {
            scal0 = sele0;
            scal1 = sele1;
         }
         osrw_du0 += scal0 * e;
         osrw_du1 += scal1 * e;


         e *= scal0;
         *eptr += e;
      }
      esum = osrw_du0;
   }
   if (vers & calc::virial) {
      for (size_t i = 0; i < virial_buffers.size(); ++i) {
         virial_buffer u = virial_buffers[i];
         virial_prec v[9];
         virial_reduce(v, u);
         double scal0 = stor0;
         double scal1 = stor1;
         if (u == vir_ev) {
            scal0 = svdw0;
            stor1 = svdw1;
         } else if (u == vir_ec || u == vir_em || u == vir_ep) {
            scal0 = sele0;
            scal1 = sele1;
         }
         for (int iv = 0; iv < 9; ++iv) {
            osrw_dv0[iv] += scal0 * v[iv];
            osrw_dv1[iv] += scal1 * v[iv];
         }
      }
      for (int iv = 0; iv < 9; ++iv) {
         vir[iv] += osrw_dv0[iv];
      }
   }
   if (vers & calc::grad) {
      double scale = stor0;
      scale_gradient(scale, gx, gy, gz);
      size_t ngrad = x_grads.size();
      for (size_t i = 1; i < ngrad; ++i) {
         grad_prec* gx1 = x_grads[i];
         grad_prec* gy1 = y_grads[i];
         grad_prec* gz1 = z_grads[i];
         double scal0 = stor0;
         double scal1 = stor1;
         if (gx1 == devx) {
            scal0 = svdw0;
            scal1 = svdw1;
         } else if (gx1 == decx || gx1 == demx || gx1 == depx) {
            scal0 = sele0;
            scal1 = sele1;
         }
         sum_gradient(scal1, osrw_dgx, osrw_dgy, osrw_dgz, gx1, gy1, gz1);
         scale_gradient(scal0, gx1, gy1, gz1);
         sum_gradient(gx, gy, gz, gx1, gy1, gz1);
      }
      sum_gradient(gx, gy, gz, osrw_gx, osrw_gy, osrw_gz);
   }
}


void osrw_energy(int vers)
{
   osrw_energy(vers, 1, default_tsconfig());
}
}
