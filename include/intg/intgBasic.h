#pragma once
#include "mdprec.h"
#include <cstdio>

namespace tinker {
class BasicIntegrator
{
protected:
   int m_iwrite;
   bool m_userat;
   void printBasic(FILE*);

public:
   bool ifSave(int istep) const;

   void rattleSave();
   void rattle(time_prec timeStep);
   void rattle2(time_prec timeStep, bool doVirial);

   BasicIntegrator();
   virtual ~BasicIntegrator();

   virtual void printDetail(FILE*) = 0;
   virtual void kickoff() = 0;
   virtual void dynamic(int istep, time_prec dt) = 0;

   void updateVelocity(time_prec t);
   void updatePosition(time_prec t);
};
}
