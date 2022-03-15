#pragma once
#include "itgbLP22.h"
#include "itgiBasic.h"

namespace tinker {
class LP22Integrator : public BasicIntegrator
{
protected:
   bool m_isNRespa1;
   const char* name() const override;
   void kickoff() override;

public:
   LP22Integrator(bool isNRespa1);
};
}
