#include "itgpLogV.h"
#include "mathfunc_sinhc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "nose.h"
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

void LogVPropagator::updatePosition(time_prec t)
{
   if (atomic) {
      if (not applyBaro) {
         propagate_pos(t);
      } else if (aniso) {
         __PlaceHolderMessage("Impl pending... updatePosition");
      } else {
         double vt = vbar * t;
         double vt2 = vt * 0.5;
         double et = std::exp(vt);
         double et2 = std::exp(vt2) * sinhc(vt2);
         propagate_pos_axbv(et, t * et2);
      }
   } else {
      __PlaceHolderMessage("Impl pending... updatePosition");
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
