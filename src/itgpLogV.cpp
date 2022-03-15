#include "itgpLogV.h"
#include "lpiston.h"
#include "mathfunc_sinhc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdppg.h"
#include "mdpq.h"
#include "nose.h"
#include "tool/trimatexp.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
//     idx  nrespa
// v1    1       1
// v2    2       1
// R0    0       1
// R1    1     nsp
// R2    2     nsp
void LogVPropagator::updateVelocityImpl(time_prec t, int idx, int nrespa)
{
   if (not applyBaro) {
      if (nrespa == 1) // v1, v2, R0
         propagate_velocity(t, gx, gy, gz);
      else // R1, R2
         propagate_velocity2(t / nrespa, gx1, gy1, gz1, t, gx2, gy2, gz2);

      return;
   }

   // if (applyBaro)
   if (atomic) {
      if (aniso) {
         double a[3][3];
         double b[3][3];
         double m[3][3];
         double tr =
            (vbar_matrix[0][0] + vbar_matrix[1][1] + vbar_matrix[2][2]) / dofP;
         for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
               m[i][j] = -vbar_matrix[i][j];
            m[i][i] -= tr;
         }

         trimat_exp(a, m, t);
         trimat_t_expm1c(b, m, t);
         // transpose
         std::swap(a[0][1], a[1][0]);
         std::swap(a[0][2], a[2][0]);
         std::swap(a[1][2], a[2][1]);
         std::swap(b[0][1], b[1][0]);
         std::swap(b[0][2], b[2][0]);
         std::swap(b[1][2], b[2][1]);

         if (nrespa == 1) // v1, v2, R0
            velAvbfAni(1, a, b, gx, gy, gz, nullptr, nullptr, nullptr);
         else // R1, R2
            velAvbfAni(nrespa, a, b, gx1, gy1, gz1, gx2, gy2, gz2);
      } else {
         double al = 1.0 + 3.0 / dofP;
         double vt = al * vbar * t;
         double vt2 = vt * 0.5;
         double a = std::exp(-vt);
         double b = t * std::exp(-vt2) * sinhc(vt2);

         if (nrespa == 1) // v1, v2, R0
            velAvbfIso(1, a, b, gx, gy, gz, nullptr, nullptr, nullptr);
         else // R1, R2
            velAvbfIso(nrespa, a, b, gx1, gy1, gz1, gx2, gy2, gz2);
      }
   } else {
      double scal[3][3], s;
      if (aniso) {
         double m[3][3];
         double tr =
            (vbar_matrix[0][0] + vbar_matrix[1][1] + vbar_matrix[2][2]) / dofP;
         for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
               m[i][j] = -vbar_matrix[i][j];
            m[i][i] -= tr;
         }
         trimat_exp(scal, m, t);
         for (int i = 0; i < 3; ++i)
            scal[i][i] -= 1;
         // transpose
         std::swap(scal[0][1], scal[1][0]);
         std::swap(scal[0][2], scal[2][0]);
         std::swap(scal[1][2], scal[2][1]);

         if (idx == 1)
            lp_propagate_mol_vel_aniso(scal);

         if (idx == 2)
            lp_propagate_mol_vel_aniso(scal);
      } else {
         double al = 1.0 + 3.0 / dofP;
         s = std::exp(-al * vbar * t) - 1;
      }

      if (idx == 1) {
         lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
         if (aniso)
            lp_propagate_mol_vel_aniso(scal);
         else
            lp_propagate_mol_vel(s);
      }

      if (nrespa == 1)
         propagate_velocity(t, gx, gy, gz);
      else
         propagate_velocity2(t / nrespa, gx1, gy1, gz1, t, gx2, gy2, gz2);

      if (idx == 2) {
         lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
         if (aniso)
            lp_propagate_mol_vel_aniso(scal);
         else
            lp_propagate_mol_vel(s);
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
         double a = std::exp(vt);
         double b = t * std::exp(vt2) * sinhc(vt2);
         propagate_pos_axbv(a, b);
      }
   } else {
      if (not applyBaro) {
         propagate_pos(t);
      } else if (aniso) {
         lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
         double scal[3][3];
         trimat_exp(scal, vbar_matrix, t);
         for (int i = 0; i < 3; ++i)
            scal[i][i] -= 1;
         pLogVPosMolAniso_acc(scal);
         propagate_pos(t);
      } else {
         lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
         double scal = std::exp(vbar * t) - 1;
         pLogVPosMolIso_acc(scal);
         propagate_pos(t);
      }
   }
}

void LogVPropagator::updateVelocity1(time_prec t)
{
   updateVelocityImpl(t, 1, 1);
}

void LogVPropagator::updateVelocity2(time_prec t)
{
   updateVelocityImpl(t, 2, 1);
}

void LogVPropagator::updateVelocityR0(time_prec t)
{
   updateVelocityImpl(t, 0, 1);
}

void LogVPropagator::updateVelocityR1(time_prec t, int nrespa)
{
   updateVelocityImpl(t, 1, nrespa);
}

void LogVPropagator::updateVelocityR2(time_prec t, int nrespa)
{
   updateVelocityImpl(t, 2, nrespa);
}
}
