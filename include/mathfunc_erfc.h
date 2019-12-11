#pragma once
#include "macro_void_cuda_def.h"
#include <cmath>


#pragma acc routine seq
__device__
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
