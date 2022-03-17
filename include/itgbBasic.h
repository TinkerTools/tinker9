#pragma once
#include "itgEnum.h"
#include "md.h"
#include <cstdio>

namespace tinker {
class BasicBarostat : public IntegratorStaticData
{
protected:
   int m_nbaro;
   void printBasic(FILE*);

public:
   BasicBarostat();
   virtual ~BasicBarostat();
   virtual void printDetail(FILE*);
   virtual BarostatEnum getBarostatEnum() const;

   virtual void control1(time_prec timeStep) {}
   virtual void control2(time_prec timeStep) {}
   virtual void control3(time_prec timeStep) {}
   virtual void control4(time_prec timeStep) {}

   virtual bool ifApply(int istep);
   bool ifApply() const;
};
}
