#include "mdpq.h"
#include "box.h"
#include "mdcalc.h"
#include "nblist.h"
#include "tool/darray.h"
#include "tool/error.h"
#include "tool/gpu_card.h"
#include <cassert>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/boxes.hh>
#include <tinker/detail/moldyn.hh>
#include <tinker/detail/usage.hh>


namespace tinker {
int rc_flag = 0;


//====================================================================//


int n;
int padded_n;
int trajn;
int nelem_buffer;


void n_data(rc_op op)
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
      nelem_buffer = gpu_max_nparallel(idevice);
      nelem_buffer = pow2_ge(nelem_buffer);
#elif TINKER_HOST
      nelem_buffer = 1;
#endif

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


void copy_pos_to_xyz(bool check_nblist)
{
   copy_pos_to_xyz_acc();
   if (check_nblist)
      refresh_neighbors();
}


void propagate_pos(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz,
                   const vel_prec* vlx, const vel_prec* vly,
                   const vel_prec* vlz)
{
   propagate_pos_acc(dt, qx, qy, qz, vlx, vly, vlz);
}


void propagate_pos(time_prec dt)
{
   propagate_pos_acc(dt, xpos, ypos, zpos, vx, vy, vz);
}


void propagate_pos_axbv(double a, double b)
{
   propagate_pos_axbv_acc(a, b);
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
         Box p;
         box_lattice(p, box_shape, l1, l2, l3, a1, a2, a3);
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

void propagate_velocity(time_prec dt, vel_prec* vlx, vel_prec* vly,
                        vel_prec* vlz, const vel_prec* vlx0,
                        const vel_prec* vly0, const vel_prec* vlz0,
                        const grad_prec* grx, const grad_prec* gry,
                        const grad_prec* grz)
{
   propagate_velocity_acc(dt, vlx, vly, vlz, vlx0, vly0, vlz0, grx, gry, grz);
}


void propagate_velocity(time_prec dt, vel_prec* vlx, vel_prec* vly,
                        vel_prec* vlz, const grad_prec* grx,
                        const grad_prec* gry, const grad_prec* grz)
{
   propagate_velocity_acc(dt, vlx, vly, vlz, grx, gry, grz);
}


void propagate_velocity(time_prec dt, const grad_prec* grx,
                        const grad_prec* gry, const grad_prec* grz)
{
   propagate_velocity(dt, vx, vy, vz, grx, gry, grz);
}


void propagate_velocity2(time_prec dt, const grad_prec* grx,
                         const grad_prec* gry, const grad_prec* grz,
                         time_prec dt2, const grad_prec* grx2,
                         const grad_prec* gry2, const grad_prec* grz2)
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


//====================================================================//


void swap_velocity(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew,
                   vel_prec* vxold, vel_prec* vyold, vel_prec* vzold)
{
   swap_velocity_acc(vxnew, vynew, vznew, vxold, vyold, vzold);
}


void propagate_pos_lp(time_prec dt, pos_prec* x_lp, pos_prec* y_lp,
                      pos_prec* z_lp, const vel_prec* vx_lp,
                      const vel_prec* vy_lp, const vel_prec* vz_lp,
                      const pos_prec* xold_lp, const pos_prec* yold_lp,
                      const pos_prec* zold_lp, double scale)
{

   propagate_pos_lp_acc(dt, x_lp, y_lp, z_lp, vx_lp, vy_lp, vz_lp, xold_lp,
                        yold_lp, zold_lp, scale);
}


void propagate_pos_lp2(time_prec dt, const pos_prec* x_lp, const pos_prec* y_lp,
                       const pos_prec* z_lp, pos_prec* xold_lp,
                       pos_prec* yold_lp, pos_prec* zold_lp, double scale)
{
   propagate_pos_lp2_acc(dt, x_lp, y_lp, z_lp, xold_lp, yold_lp, zold_lp,
                         scale);
}


void propagate_pos_lf(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz,
                      const pos_prec* qxold, const pos_prec* qyold,
                      const pos_prec* qzold, const vel_prec* vlx,
                      const vel_prec* vly, const vel_prec* vlz)
{
   propagate_pos_lf_acc(dt, qx, qy, qz, qxold, qyold, qzold, vlx, vly, vlz);
}


void propagate_velocity_lp(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
                           const vel_prec* vxnew_lp, const vel_prec* vynew_lp,
                           const vel_prec* vznew_lp, const vel_prec* vxold_lp,
                           const vel_prec* vyold_lp, const vel_prec* vzold_lp,
                           const double scale, energy_prec& eksum_new,
                           energy_prec& eksum_old)
{
   propagate_velocity_lp_acc(vx_lp, vy_lp, vz_lp, vxnew_lp, vynew_lp, vznew_lp,
                             vxold_lp, vyold_lp, vzold_lp, scale, eksum_new,
                             eksum_old);
}


void propagate_velocity_lp2(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp,
                            vel_prec* vz_lp, const pos_prec* x_lp,
                            const pos_prec* y_lp, const pos_prec* z_lp,
                            const pos_prec* xold_lp, const pos_prec* yold_lp,
                            const pos_prec* zold_lp)
{
   propagate_velocity_lp2_acc(dt, vx_lp, vy_lp, vz_lp, x_lp, y_lp, z_lp,
                              xold_lp, yold_lp, zold_lp);
}


void propagate_velocity_lp3(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
                            const vel_prec* vxnew_lp, const vel_prec* vynew_lp,
                            const vel_prec* vznew_lp, const vel_prec* vxold_lp,
                            const vel_prec* vyold_lp, const vel_prec* vzold_lp,
                            energy_prec& eksum_new)
{
   propagate_velocity_lp3_acc(vx_lp, vy_lp, vz_lp, vxnew_lp, vynew_lp, vznew_lp,
                              vxold_lp, vyold_lp, vzold_lp, eksum_new);
}
}
