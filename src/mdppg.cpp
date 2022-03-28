#include "md/pq.h"

namespace tinker {
void velAvbfIso_acc(int nrespa, vel_prec a, vel_prec b, vel_prec* vx, vel_prec* vy, vel_prec* vz,
   const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1, const grad_prec* gx2,
   const grad_prec* gy2, const grad_prec* gz2);
void mdVelAvbf(int nrespa, vel_prec a, vel_prec b, const grad_prec* gx1, const grad_prec* gy1,
   const grad_prec* gz1, const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2)
{
   velAvbfIso_acc(nrespa, a, b, vx, vy, vz, gx1, gy1, gz1, gx2, gy2, gz2);
}

void velAvbfAni_acc(int nrespa, vel_prec a[3][3], vel_prec b[3][3], vel_prec* vx, vel_prec* vy,
   vel_prec* vz, const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1,
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2);
void mdVelAvbfAn(int nrespa, vel_prec a[3][3], vel_prec b[3][3], const grad_prec* gx1,
   const grad_prec* gy1, const grad_prec* gz1, const grad_prec* gx2, const grad_prec* gy2,
   const grad_prec* gz2)
{
   velAvbfAni_acc(nrespa, a, b, vx, vy, vz, gx1, gy1, gz1, gx2, gy2, gz2);
}
}
