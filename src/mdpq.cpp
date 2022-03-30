#include "ff/box.h"
#include "ff/energy.h"
#include "ff/nblist.h"
#include "md/pq.h"
#include "tool/darray.h"
#include "tool/error.h"
#include "tool/gpucard.h"
#include <cassert>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/boxes.hh>
#include <tinker/detail/moldyn.hh>
#include <tinker/detail/usage.hh>

namespace tinker {


void mdPos(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz, const vel_prec* vlx,
   const vel_prec* vly, const vel_prec* vlz)
{
   mdPos_acc(dt, qx, qy, qz, vlx, vly, vlz);
}

void mdPos(time_prec dt)
{
   mdPos_acc(dt, xpos, ypos, zpos, vx, vy, vz);
}

void propagate_pos_axbv_acc(pos_prec a, pos_prec b);
void mdPosAxbv(pos_prec a, pos_prec b)
{
   propagate_pos_axbv_acc(a, b);
}

void mdBounds()
{
   if (!bound::use_bounds)
      return;

   mdBounds_acc();
   copyPosToXyz();
}

void mdReadFrameCopyinToXyz(std::istream& ipt, int& done)
{
   if (done)
      return;

   if (ipt) {
      std::string line;
      std::getline(ipt, line); // n and title
      std::getline(ipt, line); // either box size or first atom
      // 18.643000   18.643000   18.643000   90.000000   90.000000   90.000000
      //  1  O      8.733783    7.084710   -0.688468     1     2     3
      double l1, l2, l3, a1, a2, a3;
      int matched = std::sscanf(line.data(), "%lf%lf%lf%lf%lf%lf", &l1, &l2, &l3, &a1, &a2, &a3);
      int row = 0;
      int index;
      char name[32];
      double xr, yr, zr;
      if (matched == 6) {
         Box p;
         boxLattice(p, box_shape, l1, l2, l3, a1, a2, a3);
         boxSetCurrent(p);
      } else {
         std::sscanf(line.data(), "%d%s%lf%lf%lf", &index, name, &xr, &yr, &zr);
         index -= 1;
         atoms::x[index] = xr;
         atoms::y[index] = yr;
         atoms::z[index] = zr;
         row = 1;
      }

      for (int ir = row; ir < n; ++ir) {
         std::getline(ipt, line);
         std::sscanf(line.data(), "%d%s%lf%lf%lf", &index, name, &xr, &yr, &zr);
         index -= 1;
         atoms::x[index] = xr;
         atoms::y[index] = yr;
         atoms::z[index] = zr;
      }

      xyzData(rc_init);
   }

   if (ipt.peek() == EOF)
      done = true;
}

//====================================================================//

void mdVelB(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, //
   const vel_prec* vlx0, const vel_prec* vly0, const vel_prec* vlz0,   //
   const grad_prec* grx, const grad_prec* gry, const grad_prec* grz)
{
   mdVelB_acc(dt, vlx, vly, vlz, vlx0, vly0, vlz0, grx, gry, grz);
}

void mdVel(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz)
{
   mdVelA_acc(dt, vx, vy, vz, grx, gry, grz);
}

void mdVel2(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz,
   time_prec dt2, const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2)
{
   mdVel2_acc(dt, grx, gry, grz, dt2, grx2, gry2, grz2);
}

void mdVelData(RcOp op)
{
   if ((calc::vel & rc_flag) == 0)
      return;

   if (op & rc_dealloc) {
      darray::deallocate(vx, vy, vz);
   }

   if (op & rc_alloc) {
      darray::allocate(n, &vx, &vy, &vz);
   }

   if (op & rc_init) {
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
