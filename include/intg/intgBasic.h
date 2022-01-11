#pragma once
#include <cstdio>

namespace tinker {
class BasicIntegrator
{
protected:
   void printBasic(FILE*);

public:
   BasicIntegrator();
   virtual ~BasicIntegrator();
   virtual void printDetail(FILE*) = 0;
   virtual void dynamic(int istep, double dt) = 0;

   static void updateVelocity(double t);
   static void updatePosition(double t);
};
}
