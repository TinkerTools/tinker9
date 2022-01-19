#include "intg/intgBasic.h"
#include "mdegv.h"
#include "mdpq.h"
#include "rattle.h"
#include "tool/darray.h"
#include "tool/io_print.h"
#include <tinker/detail/inform.hh>

namespace tinker {
void BasicIntegrator::printBasic(FILE* o) {}

bool BasicIntegrator::ifSave(int istep) const
{
   return !(istep % m_iwrite);
}

void BasicIntegrator::rattleSave()
{
   if (not m_userat)
      return;

   darray::copy(g::q0, n, rattle_xold, xpos);
   darray::copy(g::q0, n, rattle_yold, ypos);
   darray::copy(g::q0, n, rattle_zold, zpos);
}

void BasicIntegrator::rattle(time_prec timeStep)
{
   if (not m_userat)
      return;

   tinker::rattle(timeStep, rattle_xold, rattle_yold, rattle_zold);
}

void BasicIntegrator::rattle2(time_prec timeStep, bool doVirial)
{
   if (not m_userat)
      return;

   tinker::rattle2(timeStep, doVirial);
}

BasicIntegrator::BasicIntegrator()
   : m_iwrite(inform::iwrite)
   , m_userat(use_rattle())
{}

BasicIntegrator::~BasicIntegrator() {}

void BasicIntegrator::updateVelocity(time_prec t)
{
   propagate_velocity(t, gx, gy, gz);
}

void BasicIntegrator::updatePosition(time_prec t)
{
   propagate_pos(t);
}
}
