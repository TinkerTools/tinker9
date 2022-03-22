#include "glob/accasync.h"
#include "math/inc.h"
#include "md/md.h"
#include "tool/darray.h"
#include "tool/error.h"

namespace tinker {
void mdDebugPosNorm_acc(pos_prec poseps, time_prec dt, //
   const vel_prec* vlx, const vel_prec* vly, const vel_prec* vlz)
{
   int which = -1;
   pos_prec tol2 = poseps * poseps;
   #pragma acc parallel loop independent async\
               copy(which) reduction(max:which)\
               deviceptr(vlx,vly,vlz)
   for (int i = 0; i < n; ++i) {
      pos_prec x1, y1, z1, norm2;
      x1 = dt * vlx[i];
      y1 = dt * vly[i];
      z1 = dt * vlz[i];
      norm2 = x1 * x1 + y1 * y1 + z1 * z1;
      bool big = norm2 > tol2;
      int flag = big ? i : -1;
      which = which > flag ? which : flag;
   }
   #pragma acc wait

   if (which >= 0) {
      printError();
      TINKER_THROW(format("MD-DEBUG POS UPDATE  --  Atom %d Tried to Move More"
                          " than %.4lf Angstroms",
         which + 1, poseps));
   }
}
}
