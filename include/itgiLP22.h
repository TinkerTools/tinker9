#pragma once
#include "itgbLP22.h"
#include "itgiBasic.h"

namespace tinker {
class LP22Integrator : public BasicIntegrator
{
protected:
   bool m_isNRespa1;
   void kickoff() override;

public:
   LP22Integrator(bool isNRespa1);
};
}
