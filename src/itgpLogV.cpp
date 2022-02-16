#include "itgpLogV.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
LogVPropagator::LogVPropagator(bool isNRespa1, bool isAtomic, bool isAniso,
                               bool isPedantic)
   : BasicPropagator()
   , m_respa(nullptr)
   , m_atomic(isAtomic)
   , m_aniso(isAniso)
   , m_pedantic(isPedantic)
{
   if (not isNRespa1) {
      nrespa = mdstuf::nrespa;
      m_respa = new RespaPropagator;
   }

   if (m_atomic)
      m_pedantic = true;
   else
      m_pedantic = false;
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
   __PlaceHolderMessage("Impl pending... updatePosition");
}

void LogVPropagator::updateVelocity0(time_prec t)
{
   __PlaceHolderMessage("Impl pending... updateVelocity0");
}

void LogVPropagator::updateVelocity1(time_prec t)
{
   __PlaceHolderMessage("Impl pending... updateVelocity1");
}

void LogVPropagator::updateVelocity2(time_prec t)
{
   __PlaceHolderMessage("Impl pending... updateVelocity2");
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
