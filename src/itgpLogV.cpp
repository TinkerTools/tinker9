#include "itgpLogV.h"
#include "lpiston.h"
#include "mathfunc_sinhc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "nose.h"
#include "tool/trimatexp.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
extern void propagate_vel_avbf_aniso_acc(double sa[3][3], double sb[3][3],
                                         const grad_prec* grx,
                                         const grad_prec* gry,
                                         const grad_prec* grz);
void LogVPropagator::updateVelocity_1_2(time_prec t, int idx)
{
   if (atomic) {
      if (applyBaro) {
         if (aniso) {
            double s1[3][3]; // v
            double s2[3][3]; // F
            double m[3][3];
            double tr =
               (vbar_matrix[0][0] + vbar_matrix[1][1] + vbar_matrix[2][2]) /
               dofP;
            for (int i = 0; i < 3; ++i) {
               for (int j = 0; j < 3; ++j)
                  m[i][j] = -vbar_matrix[i][j];
               m[i][i] -= tr;
            }

            trimat_exp(s1, m, t);
            trimat_t_expm1c(s2, m, t);
            // transpose
            std::swap(s1[0][1], s1[1][0]);
            std::swap(s1[0][2], s1[2][0]);
            std::swap(s1[1][2], s1[2][1]);
            std::swap(s2[0][1], s2[1][0]);
            std::swap(s2[0][2], s2[2][0]);
            std::swap(s2[1][2], s2[2][1]);

            propagate_vel_avbf_aniso_acc(s1, s2, gx, gy, gz);
         } else {
            double al = 1.0 + 3.0 / dofP;
            double vt = al * vbar * t;
            double vt2 = vt * 0.5;
            double et = std::exp(-vt);
            double et2 = std::exp(-vt2) * sinhc(vt2);
            propagate_vel_avbf(et, t * et2, gx, gy, gz);
         }
      } else {
         if (idx == 1)
            BasicPropagator::updateVelocity1(t);
         else if (idx == 2)
            BasicPropagator::updateVelocity2(t);
      }
   } else {
      if (applyBaro) {
         if (aniso) {
            double scal[3][3];
            double m[3][3];
            double tr =
               (vbar_matrix[0][0] + vbar_matrix[1][1] + vbar_matrix[2][2]) /
               dofP;
            for (int i = 0; i < 3; ++i) {
               for (int j = 0; j < 3; ++j)
                  m[i][j] = vbar_matrix[i][j];
               m[i][i] += tr;
            }
            trimat_exp(scal, m, -t);
            // transpose
            std::swap(scal[0][1], scal[1][0]);
            std::swap(scal[0][2], scal[2][0]);
            std::swap(scal[1][2], scal[2][1]);

            if (idx == 1) {
               lp_propagate_mol_vel_aniso(scal);
               BasicPropagator::updateVelocity1(t);
            } else if (idx == 2) {
               BasicPropagator::updateVelocity2(t);
               lp_propagate_mol_vel_aniso(scal);
            }
         } else {
            double al = 1.0 + 3.0 / dofP;
            double scal = std::exp(al * vbar * t) - 1;

            if (idx == 1) {
               lp_propagate_mol_vel(scal);
               BasicPropagator::updateVelocity1(t);
            } else if (idx == 2) {
               BasicPropagator::updateVelocity2(t);
               lp_propagate_mol_vel(scal);
            }
         }
      } else {
         if (idx == 1)
            BasicPropagator::updateVelocity1(t);
         else if (idx == 2)
            BasicPropagator::updateVelocity2(t);
      }
   }
}

LogVPropagator::LogVPropagator(bool isNRespa1)
   : BasicPropagator()
   , m_respa(nullptr)
{
   if (not isNRespa1) {
      nrespa = mdstuf::nrespa;
      m_respa = new RespaPropagator;
   }
}

LogVPropagator::~LogVPropagator()
{
   if (m_respa) {
      delete m_respa;
      m_respa = nullptr;
   }
}

extern void pLogVPosMolIso_acc(double scal);
extern void pLogVPosMolAniso_acc(double (*scal)[3]);
extern void pLogVPosAtmAniso_acc(double (*a)[3], double (*b)[3]);
void LogVPropagator::updatePosition(time_prec t)
{
   if (atomic) {
      if (not applyBaro) {
         propagate_pos(t);
      } else if (aniso) {
         double a[3][3], b[3][3];
         trimat_exp(a, vbar_matrix, t);
         trimat_t_expm1c(b, vbar_matrix, t);
         pLogVPosAtmAniso_acc(a, b);
      } else {
         double vt = vbar * t;
         double vt2 = vt * 0.5;
         double et = std::exp(vt);
         double et2 = std::exp(vt2) * sinhc(vt2);
         propagate_pos_axbv(et, t * et2);
      }
   } else {
      if (not applyBaro) {
         propagate_pos(t);
      } else if (aniso) {
         lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
         double scal[3][3];
         trimat_exp(scal, vbar_matrix, t);
         pLogVPosMolAniso_acc(scal);
         propagate_pos(t);
      } else {
         lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
         double scal = std::exp(vbar * t);
         pLogVPosMolIso_acc(scal);
         propagate_pos(t);
      }
   }
}

void LogVPropagator::updateVelocity0(time_prec t)
{
   __PlaceHolderMessage("Impl pending... updateVelocity0");
}

void LogVPropagator::updateVelocity1(time_prec t)
{
   updateVelocity_1_2(t, 1);
}

void LogVPropagator::updateVelocity2(time_prec t)
{
   updateVelocity_1_2(t, 2);
}

void LogVPropagator::updateVelocityR1(time_prec tfast, time_prec t)
{
   __PlaceHolderMessage("Impl pending... updateVelocityR1");
}

void LogVPropagator::updateVelocityR2(time_prec tfast, time_prec t)
{
   __PlaceHolderMessage("Impl pending... updateVelocityR2");
}
}
