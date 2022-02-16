#pragma once
#include "itgiBasic.h"

namespace tinker {
class Nhc06Integrator : public BasicIntegrator
{
protected:
   bool m_isNRespa1;
   void kickoff() override;

public:
   Nhc06Integrator(bool isNRespa1);
};
}
