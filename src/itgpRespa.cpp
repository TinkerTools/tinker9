#include "itgpRespa.h"
#include "md.h"
#include "tool/darray.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
RespaPropagator::RespaPropagator()
{
   darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);
   nrespa = mdstuf::nrespa;
}

RespaPropagator::~RespaPropagator()
{
   darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);
}

void RespaPropagator::updateVelocityR1(time_prec t, int nrespa)
{
   mdVel2(t / nrespa, gx1, gy1, gz1, t, gx2, gy2, gz2);
}

void RespaPropagator::updateVelocityR2(time_prec t, int nrespa)
{
   this->updateVelocityR1(t, nrespa);
}
}
