#include "ff/energy.h"
#include "md/integrator.h"
#include "md/misc.h"
#include "md/pq.h"
#include "tool/darray.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
RespaDevice::RespaDevice()
{
   darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);
   nrespa = mdstuf::nrespa;
}

RespaDevice::~RespaDevice()
{
   darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);
}

void RespaDevice::velR0(time_prec t)
{
   mdVel(t, gx, gy, gz);
}

void RespaDevice::velR1(time_prec t, int nrespa)
{
   mdVel2(t / nrespa, gx1, gy1, gz1, t, gx2, gy2, gz2);
}

void RespaDevice::velR2(time_prec t, int nrespa)
{
   this->velR1(t, nrespa);
}
}
