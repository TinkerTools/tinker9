
#include "mdpq.h"
#include "box.h"
#include "darray.h"
#include "error.h"
#include "gpu_card.h"
#include "mdcalc.h"
#include "nblist.h"
#include "tinker_rt.h"
#include <cassert>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/boxes.hh>
#include <tinker/detail/moldyn.hh>
#include <tinker/detail/usage.hh>


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


      if (usage::nuse != n) {
         TINKER_THROW("All atoms must be active.");
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
      refresh_neighbors();
}


void bounds()
{
   if (!bound::use_bounds)
      return;


   bounds_pos_acc();
   copy_pos_to_xyz();
}


void read_frame_copyin_to_xyz(std::istream& ipt, int& done)
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
      int matched = std::sscanf(line.data(), "%lf%lf%lf%lf%lf%lf", &l1, &l2,
                                &l3, &a1, &a2, &a3);
      int row = 0;
      int index;
      char name[32];
      double xr, yr, zr;
      if (matched == 6) {
         boxes::xbox = l1;
         boxes::ybox = l2;
         boxes::zbox = l3;
         boxes::alpha = a1;
         boxes::beta = a2;
         boxes::gamma = a3;
         TINKER_RT(lattice)();

         Box p;
         get_tinker_box_module(p);
         set_default_box(p);
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

      xyz_data(rc_init);
   }


   if (ipt.peek() == EOF)
      done = true;
}


void xyz_data(rc_op op)
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
         darray::copyin(PROCEED_NEW_Q, n, x, atoms::x);
         darray::copyin(PROCEED_NEW_Q, n, y, atoms::y);
         darray::copyin(PROCEED_NEW_Q, n, z, atoms::z);
      } else {
         darray::copyin(PROCEED_NEW_Q, n, xpos, atoms::x);
         darray::copyin(PROCEED_NEW_Q, n, ypos, atoms::y);
         darray::copyin(PROCEED_NEW_Q, n, zpos, atoms::z);
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
      darray::deallocate(mass, massinv);
   }

   if (op & rc_alloc) {
      darray::allocate(n, &mass, &massinv);
   }

   if (op & rc_init) {
      std::vector<double> mbuf(n);
      for (int i = 0; i < n; ++i)
         mbuf[i] = 1 / atomid::mass[i];
      darray::copyin(PROCEED_NEW_Q, n, massinv, mbuf.data());
      darray::copyin(WAIT_NEW_Q, n, mass, atomid::mass);
   }
}


void vel_data(rc_op op)
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
      darray::copyin2(PROCEED_NEW_Q, 0, 3, n, vx, moldyn::v);
      darray::copyin2(PROCEED_NEW_Q, 1, 3, n, vy, moldyn::v);
      darray::copyin2(WAIT_NEW_Q, 2, 3, n, vz, moldyn::v);
   }
}
TINKER_NAMESPACE_END
