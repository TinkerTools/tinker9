#pragma once
#include "mdprec.h"
#include <cstdio>

namespace tinker {
enum class BarostatEnum;

class BasicBarostat
{
protected:
   BarostatEnum m_baroEnum;
   int m_nbaro;
   void printBasic(FILE*);
   BasicBarostat(BarostatEnum);

public:
   BasicBarostat();
   virtual ~BasicBarostat();
   virtual void printDetail(FILE*);
   virtual void control1(time_prec timeStep) {}
   virtual void control2(time_prec timeStep) {}
   virtual void control3(time_prec timeStep) {}
   virtual void control4(time_prec timeStep) {}
   virtual bool ifApply(int istep);
   BarostatEnum getBarostatEnum() const;
};
}
