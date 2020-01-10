#pragma once
#include "macro.h"
#include "seq_def.h"
#include <cmath>


TINKER_NAMESPACE_BEGIN
SEQ_ROUTINE
inline float erfcf_hastings(float x)
{
   float exp2a = expf(-x * x);
   float t = 1.0f / (1.0f + 0.3275911f * x);
   return (0.254829592f +
           (-0.284496736f +
            (1.421413741f + (-1.453152027f + 1.061405429f * t) * t) * t) *
              t) *
      t * exp2a;
}
TINKER_NAMESPACE_END
