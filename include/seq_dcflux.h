#pragma once
#include "mathfunc.h"
#include "realn.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
SEQ_CUDA
inline void bnd_dcflux(real* restrict dcfa, real* restrict dcfb, real ra,
                       real rb, real pb, real pota, real potb)
{
   real3 dba = ra - rb;
   real rab = REAL_SQRT(dba.x * dba.x + dba.y * dba.y + dba.z * dba.z);
   real dpot = pota - potb;
   pb = pb / rab;

   real3 ddqr = pb * dr;
   dcfa += dpot * ddqr;
   dcfb -= dpot * ddqr;
   
}

#pragma acc routine seq
SEQ_CUDA
inline void ang_dcflux(real* restrict dcfa, real* restrict dcfb,
                       real* restrict dcfc, real ra, real rb, real rc, real pa1,
                       real pa2, real pb1, real pb2, real pota, real potb,
                       real potc)
{
   real3 dba = ra - rb;
   real3 dbc = rc - rb;
   real rba2 = dba.x * dba.x + dba.y * dba.y + dba.z * dba.z;
   real rba = REAL_SQRT(rba2);
   real rba3 = rba2 * rba;

   real rbc2 = dbc.x * dbc.x + dbc.y * dbc.y + dbc.z * dbc.z;
   real rbc = REAL_SQRT(rbc2);
   real rbc3 = rbc2 * rbc;

   real dpota = pota - potb;
   real dpotc = potc - potb;
   pb1 *= dpota;
   pb2 *= dpotc;

   real c1 = pb2 / rba;
   real c2 = pb1 / rbc;
   real3 fa1 = c1 * dba;
   real3 fc1 = c2 * dbc;
   real3 fb1 = -fa1 - fc1;

   dcfa += fa1;
   dcfb += fb1;
   dcfc += fc1;

   // real fxa1 = pb2 * dba.x / rba;
   // real fya1 = pb2 * dba.y / rba;
   // real fza1 = pb2 * dba.z / rba;
   // real fxc1 = pb1 * dbc.x / rbc;
   // real fyc1 = pb1 * dbc.y / rbc;
   // real fzc1 = pb1 * dbc.z / rbc;
   // real fxb1 = -fxa1 - fxc1;
   // real fyb1 = -fya1 - fyc1;
   // real fzb1 = -fza1 - fzc1;

   real dot = dot3(dba, dbc);
   real term = -rba * rbc / REAL_SQRT(rba2 * rbc2 - dot * dot);
   real fterm = term * (dpota * pa1 + dpotc * pa2);

   c1 = 1 / (rba * rbc);
   c2 = dot / (rba3 * rbc);
   real c3 = dot / (rbc3 * rba);
   real3 terma = c1 * dbc - c2 * dba;
   real3 termc = c1 * dba - c3 * dbc;

   real3 fa2 = fterm * terma;
   real3 fc2 = fterm * termc;
   real3 fb2 = -fa2 - fc2;

   dcfa += fa2;
   dcfb += fb2;
   dcfc += fc2;
}
}
