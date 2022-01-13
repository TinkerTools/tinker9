#pragma once
#include "mdprec.h"
#include <cstdio>

namespace tinker {
class BasicBarostat
{
protected:
   int m_nbaro;
   void printBasic(FILE*);

public:
   BasicBarostat();
   virtual ~BasicBarostat();
   virtual void printDetail(FILE*);
   virtual void control1(time_prec timeStep) {}
   virtual void control2(time_prec timeStep) {}
   virtual void control3(time_prec timeStep) {}
   virtual void control4(time_prec timeStep) {}
   virtual bool ifApply(int istep);
};
}
