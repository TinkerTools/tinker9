
#include "md_pq.h"
#include "dev_array.h"
#include "gpu_card.h"
#include "md_calc.h"
#include "nblist.h"
#include <cassert>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/moldyn.hh>


TINKER_NAMESPACE_BEGIN
int rc_flag = 0;


//====================================================================//


int n;
int padded_n;
int trajn;


void n_data(rc_op op)
{
   if (op & rc_dealloc) {
      trajn = -1;
      n = 0;
      padded_n = 0;
   }

   if (op & rc_alloc) {
      n = atoms::n;
      padded_n = (n + WARP_SIZE - 1) / WARP_SIZE;
      padded_n *= WARP_SIZE;

      if (calc::traj & rc_flag) {
         // trajn must have been initialized by this point
         assert(trajn >= 0);
      }
   }
}


//====================================================================//


real *x, *y, *z;
real *trajx, *trajy, *trajz;
pos_prec *xpos, *ypos, *zpos;


void copy_pos_to_xyz()
{
   copy_pos_to_xyz_acc();
}


void propagate_xyz(time_prec dt, bool check_nblist)
{
   propagate_pos_acc(dt);
   copy_pos_to_xyz();
   if (check_nblist)
      nblist_data(rc_evolve);
}


void xyz_data(rc_op op)
{
   if ((calc::xyz & rc_flag) == 0)
      return;

   if (op & rc_dealloc) {
      if (calc::traj & rc_flag) {
         device_array::deallocate(trajx, trajy, trajz);
         x = nullptr;
         y = nullptr;
         z = nullptr;
      } else {
         trajx = nullptr;
         trajy = nullptr;
         trajz = nullptr;
         device_array::deallocate(x, y, z);
         if (sizeof(pos_prec) == sizeof(real)) {
            xpos = nullptr;
            ypos = nullptr;
            zpos = nullptr;
         } else {
            device_array::deallocate(xpos, ypos, zpos);
         }
      }
   }

   if (op & rc_alloc) {
      if (calc::traj & rc_flag) {
         device_array::allocate(n * trajn, &trajx, &trajy, &trajz);
         x = trajx;
         y = trajy;
         z = trajz;
      } else {
         device_array::allocate(n, &x, &y, &z);
         if (sizeof(pos_prec) == sizeof(real)) {
            xpos = (pos_prec*)x;
            ypos = (pos_prec*)y;
            zpos = (pos_prec*)z;
         } else {
            device_array::allocate(n, &xpos, &ypos, &zpos);
         }
      }
   }

   if (op & rc_init) {
      if (calc::traj & rc_flag) {
         device_array::copyin(PROCEED_NEW_Q, n, x, atoms::x);
         device_array::copyin(PROCEED_NEW_Q, n, y, atoms::y);
         device_array::copyin(PROCEED_NEW_Q, n, z, atoms::z);
      } else {
         device_array::copyin(PROCEED_NEW_Q, n, xpos, atoms::x);
         device_array::copyin(PROCEED_NEW_Q, n, ypos, atoms::y);
         device_array::copyin(PROCEED_NEW_Q, n, zpos, atoms::z);
         copy_pos_to_xyz();
      }
   }
}


//====================================================================//


mass_prec *mass, *massinv;


vel_prec *vx, *vy, *vz;


void propagate_velocity(time_prec dt, const real* grx, const real* gry,
                        const real* grz)
{
   propagate_velocity_acc(dt, grx, gry, grz);
}


void propagate_velocity(time_prec dt, const fixed* grx, const fixed* gry,
                        const fixed* grz)
{
   propagate_velocity_acc(dt, grx, gry, grz);
}


void propagate_velocity2(time_prec dt, const real* grx, const real* gry,
                         const real* grz, time_prec dt2, const real* grx2,
                         const real* gry2, const real* grz2)
{
   propagate_velocity2_acc(dt, grx, gry, grz, dt2, grx2, gry2, grz2);
}


void propagate_velocity2(time_prec dt, const fixed* grx, const fixed* gry,
                         const fixed* grz, time_prec dt2, const fixed* grx2,
                         const fixed* gry2, const fixed* grz2)
{
   propagate_velocity2_acc(dt, grx, gry, grz, dt2, grx2, gry2, grz2);
}


void mass_data(rc_op op)
{
   if ((calc::mass & rc_flag) == 0)
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(mass, massinv);
   }

   if (op & rc_alloc) {
      device_array::allocate(n, &mass, &massinv);
   }

   if (op & rc_init) {
      std::vector<double> mbuf(n);
      for (int i = 0; i < n; ++i)
         mbuf[i] = 1 / atomid::mass[i];
      device_array::copyin(PROCEED_NEW_Q, n, massinv, mbuf.data());
      device_array::copyin(WAIT_NEW_Q, n, mass, atomid::mass);
   }
}


void vel_data(rc_op op)
{
   if ((calc::vel & rc_flag) == 0)
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(vx, vy, vz);
   }

   if (op & rc_alloc) {
      device_array::allocate(n, &vx, &vy, &vz);
   }

   if (op & rc_init) {
      device_array::copyin2(PROCEED_NEW_Q, 0, 3, n, vx, moldyn::v);
      device_array::copyin2(PROCEED_NEW_Q, 1, 3, n, vy, moldyn::v);
      device_array::copyin2(WAIT_NEW_Q, 2, 3, n, vz, moldyn::v);
   }
}
TINKER_NAMESPACE_END
