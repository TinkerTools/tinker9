#pragma once
#include "ff/precision.h"

namespace tinker {
// udir = polarity * field

__global__
void pcgUdirV1(int n, const real* polarity, //
   real (*udir)[3], const real (*field)[3]);

__global__
void pcgUdirV2(int n, const real* polarity, //
   real (*udir)[3], real (*udirp)[3], const real (*field)[3], const real (*fieldp)[3]);

// r(0) = E - (1/polarity + Tu) u(0) = (udir - u(0))/polarity + mutual field

__global__
void pcgRsd0V1(int n, const real* polarity_inv, real (*rsd)[3], //
   const real (*udir)[3], const real (*uind)[3], const real (*field)[3]);

__global__
void pcgRsd0V2(int n, const real* polarity_inv, real (*rsd)[3], real (*rsp)[3], //
   const real (*udir)[3], const real (*udip)[3], const real (*uind)[3], const real (*uinp)[3],
   const real (*field)[3], const real (*fielp)[3]);

__global__
void pcgRsd0V3(int n, const real* polarity_inv, real (*rsd)[3], //
   const real (*udir)[3], const real (*uind)[3], const real (*field)[3],
   const real (*polscale)[3][3]);
}
