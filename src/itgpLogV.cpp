#include "itgpLogV.h"
#include "lpiston.h"
#include "mathfunc_sinhc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "nose.h"
#include "tool/trimatexp.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
void LogVPropagator::updateVelocity_1_2(time_prec t, int idx)
{
   if (atomic) {
      if (applyBaro) {
         if (aniso) {
            __PlaceHolderMessage("Impl pending... updateVelocity1_2");
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
         else
            BasicPropagator::updateVelocity2(t);
      }
   } else {
      __PlaceHolderMessage("Impl pending... updateVelocity1_2");
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
