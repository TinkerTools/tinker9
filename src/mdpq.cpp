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
void nData(RcOp op)
{
   if (op & rc_dealloc) {
      trajn = -1;
      n = 0;
      padded_n = 0;
      nelem_buffer = 0;
   }

   if (op & rc_alloc) {
      n = atoms::n;
      padded_n = (n + WARP_SIZE - 1) / WARP_SIZE;
      padded_n *= WARP_SIZE;

      if (calc::traj & rc_flag) {
         // trajn must have been initialized by this point
         assert(trajn >= 0);
      }

#if TINKER_CUDART
      nelem_buffer = gpuMaxNParallel(idevice);
      nelem_buffer = pow2Ge(nelem_buffer);
#elif TINKER_HOST
      nelem_buffer = 1;
#endif

      if (usage::nuse != n) {
         TINKER_THROW("All atoms must be active.");
      }
   }
}

//====================================================================//

void mdCopyPosToXyz()
{
   mdCopyPosToXyz_acc();
}

void mdCopyPosToXyz(bool check_nblist)
{
   mdCopyPosToXyz_acc();
   if (check_nblist)
      nblistRefresh();
}

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
   mdCopyPosToXyz();
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

void xyzData(RcOp op)
{
   if ((calc::xyz & rc_flag) == 0)
      return;

   if (op & rc_dealloc) {
      if (calc::traj & rc_flag) {
         darray::deallocate(trajx, trajy, trajz);
         x = nullptr;
         y = nullptr;
         z = nullptr;
      } else {
         trajx = nullptr;
         trajy = nullptr;
         trajz = nullptr;
         darray::deallocate(x, y, z);
         if (sizeof(pos_prec) == sizeof(real)) {
            xpos = nullptr;
            ypos = nullptr;
            zpos = nullptr;
         } else {
            darray::deallocate(xpos, ypos, zpos);
         }
      }
   }

   if (op & rc_alloc) {
      if (calc::traj & rc_flag) {
         darray::allocate(n * trajn, &trajx, &trajy, &trajz);
         x = trajx;
         y = trajy;
         z = trajz;
      } else {
         darray::allocate(n, &x, &y, &z);
         if (sizeof(pos_prec) == sizeof(real)) {
            xpos = (pos_prec*)x;
            ypos = (pos_prec*)y;
            zpos = (pos_prec*)z;
         } else {
            darray::allocate(n, &xpos, &ypos, &zpos);
         }
      }
   }

   if (op & rc_init) {
      if (calc::traj & rc_flag) {
         darray::copyin(g::q0, n, x, atoms::x);
         darray::copyin(g::q0, n, y, atoms::y);
         darray::copyin(g::q0, n, z, atoms::z);
      } else {
         darray::copyin(g::q0, n, xpos, atoms::x);
         darray::copyin(g::q0, n, ypos, atoms::y);
         darray::copyin(g::q0, n, zpos, atoms::z);
         mdCopyPosToXyz();
      }
   }
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

void massData(RcOp op)
{
   if ((calc::mass & rc_flag) == 0)
      return;

   if (op & rc_dealloc) {
      darray::deallocate(mass, massinv);
   }

   if (op & rc_alloc) {
      darray::allocate(n, &mass, &massinv);
   }

   if (op & rc_init) {
      std::vector<double> mbuf(n);
      for (int i = 0; i < n; ++i)
         mbuf[i] = 1 / atomid::mass[i];
      darray::copyin(g::q0, n, massinv, mbuf.data());
      darray::copyin(g::q0, n, mass, atomid::mass);
      wait_for(g::q0);
   }
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
      wait_for(g::q0);
   }
}
}
