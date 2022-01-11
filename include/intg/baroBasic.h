#pragma once
#include <cstdio>

namespace tinker {
class BasicBarostat
{
protected:
   int m_nbaro;
   void printBasic(FILE*);

public:
   int nbaro() const;
   BasicBarostat();
   virtual ~BasicBarostat();
   virtual void printDetail(FILE*);
   virtual void control(double timeStep);
   virtual void controlAfter(double timeStep);
};
}
