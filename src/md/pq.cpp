#include "md/pq.h"
#include "tool/darray.h"
#include <tinker/detail/moldyn.hh>

namespace tinker {
extern void mdPos_acc(time_prec, pos_prec*, pos_prec*, pos_prec*, //
   const vel_prec*, const vel_prec*, const vel_prec*);
void mdPos(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz, const vel_prec* vlx,
   const vel_prec* vly, const vel_prec* vlz)
{
   mdPos_acc(dt, qx, qy, qz, vlx, vly, vlz);
}

void mdPos(time_prec dt)
{
   mdPos_acc(dt, xpos, ypos, zpos, vx, vy, vz);
}

extern void mdPosAxbv_acc(pos_prec a, pos_prec b);
void mdPosAxbv(pos_prec a, pos_prec b)
{
   mdPosAxbv_acc(a, b);
}

extern void mdPosAxbvAn_acc(pos_prec (*a)[3], pos_prec (*b)[3]);
void mdPosAxbvAn(pos_prec (*a)[3], pos_prec (*b)[3])
{
   mdPosAxbvAn_acc(a, b);
}

extern void mdDebugPosNorm_acc(pos_prec poseps, time_prec dt, //
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void mdDebugPosNorm(pos_prec poseps, time_prec dt, //
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
   mdDebugPosNorm_acc(poseps, dt, vx, vy, vz);
}
}

namespace tinker {
extern void mdVel_acc(time_prec, vel_prec*, vel_prec*, vel_prec*, //
   const grad_prec*, const grad_prec*, const grad_prec*);
void mdVel(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz)
{
   mdVel_acc(dt, vx, vy, vz, grx, gry, grz);
}

extern void mdVelB_acc(time_prec, vel_prec*, vel_prec*, vel_prec*, //
   const vel_prec*, const vel_prec*, const vel_prec*,              //
   const grad_prec*, const grad_prec*, const grad_prec*);
void mdVelB(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, //
   const vel_prec* vlx0, const vel_prec* vly0, const vel_prec* vlz0,   //
   const grad_prec* grx, const grad_prec* gry, const grad_prec* grz)
{
   mdVelB_acc(dt, vlx, vly, vlz, vlx0, vly0, vlz0, grx, gry, grz);
}

extern void mdVel2_acc(time_prec, const grad_prec*, const grad_prec*, const grad_prec*, time_prec,
   const grad_prec*, const grad_prec*, const grad_prec*);
void mdVel2(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz,
   time_prec dt2, const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2)
{
   mdVel2_acc(dt, grx, gry, grz, dt2, grx2, gry2, grz2);
}

extern void mdVelAvbf_acc(int nrespa, vel_prec a, vel_prec b, vel_prec* vx, vel_prec* vy,
   vel_prec* vz, const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1,
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2);
void mdVelAvbf(int nrespa, vel_prec a, vel_prec b, const grad_prec* gx1, const grad_prec* gy1,
   const grad_prec* gz1, const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2)
{
   mdVelAvbf_acc(nrespa, a, b, vx, vy, vz, gx1, gy1, gz1, gx2, gy2, gz2);
}

extern void mdVelAvbfAn_acc(int nrespa, vel_prec a[3][3], vel_prec b[3][3], vel_prec* vx,
   vel_prec* vy, vel_prec* vz, const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1,
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2);
void mdVelAvbfAn(int nrespa, vel_prec a[3][3], vel_prec b[3][3], const grad_prec* gx1,
   const grad_prec* gy1, const grad_prec* gz1, const grad_prec* gx2, const grad_prec* gy2,
   const grad_prec* gz2)
{
   mdVelAvbfAn_acc(nrespa, a, b, vx, vy, vz, gx1, gy1, gz1, gx2, gy2, gz2);
}
}

namespace tinker {
void mdVelData(RcOp op)
{
   if (not(calc::vel & rc_flag))
      return;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(vx, vy, vz);
   }

   if (op & RcOp::ALLOC) {
      darray::allocate(n, &vx, &vy, &vz);
   }

   if (op & RcOp::INIT) {
      std::vector<vel_prec> vvx(n), vvy(n), vvz(n);
      for (int i = 0; i < n; ++i) {
         vvx[i] = moldyn::v[3 * i + 0];
         vvy[i] = moldyn::v[3 * i + 1];
         vvz[i] = moldyn::v[3 * i + 2];
      }
      darray::copyin(g::q0, n, vx, vvx.data());
      darray::copyin(g::q0, n, vy, vvy.data());
      darray::copyin(g::q0, n, vz, vvz.data());
      waitFor(g::q0);
   }
}
}
