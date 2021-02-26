#pragma once
#include "glob.nelembuffer.h"
#include "mdprec.h"
#include "tool/rc_man.h"
#include <istream>


namespace tinker {
/**
 * \ingroup mdcalc
 * \var rc_flag
 * \brief Global bitmask.
 */
extern int rc_flag;


//====================================================================//


/**
 * \ingroup mdpq
 * \{
 * \var n
 * \brief Number of atoms.
 *
 * \var padded_n
 * \brief Number of atoms padded by #WARP_SIZE.
 * \see WARP_SIZE
 *
 * \var trajn
 * \brief Number of the trajectory frames.
 * \}
 */
extern int n;
extern int padded_n;
extern int trajn;


void n_data(rc_op);


//====================================================================//


/**
 * \ingroup mdpq
 * \{
 * \var x
 * \brief Current coordinates used in energy evaluation and neighbor lists.
 * \var y
 * \copydoc x
 * \var z
 * \copydoc x
 *
 * \var trajx
 * \brief Coordinates of all the trajectory frames.
 * \var trajy
 * \copydoc trajx
 * \var trajz
 * \copydoc trajx
 *
 * \var xpos
 * \brief Coordinates used in integrators.
 * \note
 *    - New arrays will be allocated only if `sizeof(pos_prec) > sizeof(real)`,
 *    otherwise, they will be aliases of #x, #y, and #z.
 *    - Whenever #xpos, #ypos, #zpos get updated by integrators, barostats etc.,
 *    #x, #y, #z must be updated immediately.
 * \see pos_prec
 * \see real
 * \var ypos
 * \copydoc xpos
 * \var zpos
 * \copydoc xpos
 * \}
 */
extern real *x, *y, *z;
extern real *trajx, *trajy, *trajz;
extern pos_prec *xpos, *ypos, *zpos;
static_assert(sizeof(pos_prec) >= sizeof(real),
              "Type pos_prec cannot be shorter than type real.");


/**
 * \ingroup mdpq
 * \brief Update #x, #y, #z by #xpos, #ypos, and #zpos.
 * If #xpos etc. are only aliases, return directly.
 */
void copy_pos_to_xyz();
void copy_pos_to_xyz(bool check_nblist);
void copy_pos_to_xyz_acc();


/**
 * \ingroup mdpq
 * \brief Updates #x, #y, #z via `x += v * dt`.
 * Currently #xpos, #ypos, and #zpos are integrated first, then uses
 * #copy_pos_to_xyz() to update #x, #y, #z.
 * \param dt            Time-step for this update.
 * \param check_nblist  If `ture`, check the neighbor lists after updating the
 *                      coordinates.
 */
void propagate_pos(time_prec, pos_prec*, pos_prec*, pos_prec*, const vel_prec*,
                   const vel_prec*, const vel_prec*);
void propagate_pos(time_prec);
void propagate_pos_acc(time_prec, pos_prec*, pos_prec*, pos_prec*,
                       const vel_prec*, const vel_prec*, const vel_prec*);


/**
 * \ingroup mdpq
 * x = a x + b v
 */
void propagate_pos_axbv(double a, double b);
void propagate_pos_axbv_acc(double a, double b);


/**
 * \ingroup mdpq
 * \brief Call #bounds() at least every x steps in MD.
 */
constexpr int BOUNDS_EVERY_X_STEPS = 500;
/**
 * \ingroup mdpq
 * \brief Finds the geometric center of each molecule and translate any stray
 * molecules back into the periodic box on GPU.
 * \note
 *    - Updating #x, #y, #z is the goal.
 *    - Checks whether PBC is in use inside this function.
 *    - Will not perturb the neighbor lists so no need to update them.
 *    - Tinker uses centers of mass.
 */
void bounds();
void bounds_pos_acc();


void read_frame_copyin_to_xyz(std::istream& input, int& done);
void xyz_data(rc_op);


//====================================================================//


/**
 * \ingroup mdpq
 * \{
 * \var mass
 * \brief Atomic mass.
 * \var massinv
 * \brief Inversed atomic mass.
 *
 * \var vx
 * \brief Velocities.
 * \var vy
 * \copydoc vx
 * \var vz
 * \copydoc vx
 * \}
 */
extern mass_prec* mass;
extern mass_prec* massinv;
extern vel_prec *vx, *vy, *vz;


// v = v0 + g/m dt; v and v0 are different
void propagate_velocity(time_prec dt, vel_prec* vlx, vel_prec* vly,
                        vel_prec* vlz, const vel_prec* vlx0,
                        const vel_prec* vly0, const vel_prec* vlz0,
                        const grad_prec* grx, const grad_prec* gry,
                        const grad_prec* grz);
void propagate_velocity_acc(time_prec, vel_prec*, vel_prec*, vel_prec*,
                            const vel_prec*, const vel_prec*, const vel_prec*,
                            const grad_prec*, const grad_prec*,
                            const grad_prec*);


void propagate_velocity(time_prec dt, vel_prec* vlx, vel_prec* vly,
                        vel_prec* vlz, const grad_prec* grx,
                        const grad_prec* gry, const grad_prec* grz);
void propagate_velocity_acc(time_prec, vel_prec*, vel_prec*, vel_prec*,
                            const grad_prec*, const grad_prec*,
                            const grad_prec*);


/**
 * \ingroup mdpq
 * \brief Update velocities via `v += -g/m dt`.
 */
void propagate_velocity(time_prec dt, const grad_prec* grx,
                        const grad_prec* gry, const grad_prec* grz);
/**
 * \ingroup mdpq
 * \brief Update velocities via `v += (-g/m dt -g2/m dt2)`.
 */
void propagate_velocity2(time_prec dt, const grad_prec* grx,
                         const grad_prec* gry, const grad_prec* grz,
                         time_prec dt2, const grad_prec* grx2,
                         const grad_prec* gry2, const grad_prec* grz2);
void propagate_velocity2_acc(time_prec, const grad_prec*, const grad_prec*,
                             const grad_prec*, time_prec, const grad_prec*,
                             const grad_prec*, const grad_prec*);


void mass_data(rc_op);
void vel_data(rc_op);


//====================================================================//


void swap_velocity(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew,
                   vel_prec* vxold, vel_prec* vyold, vel_prec* vzold);
void swap_velocity_acc(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew,
                       vel_prec* vxold, vel_prec* vyold, vel_prec* vzold);


void propagate_pos_lp(time_prec dt, pos_prec* x_lp, pos_prec* y_lp,
                      pos_prec* z_lp, const vel_prec* vx_lp, const vel_prec* vy,
                      const vel_prec* vz, const pos_prec* xold_lp,
                      const pos_prec* yold_lp, const pos_prec* zold_lp,
                      double scale);
void propagate_pos_lp_acc(time_prec dt, pos_prec* x_lp, pos_prec* y_lp,
                          pos_prec* z_lp, const vel_prec* vx_lp,
                          const vel_prec* vy_lp, const vel_prec* vz_lp,
                          const pos_prec* xold_lp, const pos_prec* yold_lp,
                          const pos_prec* zold_lp, double scale);


void propagate_pos_lp2(time_prec dt, const pos_prec* x_lp, const pos_prec* y_lp,
                       const pos_prec* z_lp, pos_prec* xold_lp,
                       pos_prec* yold_lp, pos_prec* zold_lp, double scale);
void propagate_pos_lp2_acc(time_prec dt, const pos_prec* x_lp,
                           const pos_prec* y_lp, const pos_prec* z_lp,
                           pos_prec* xold_lp, pos_prec* yold_lp,
                           pos_prec* zold_lp, double scale);


void propagate_pos_lf(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz,
                      const pos_prec* qxold, const pos_prec* qyold,
                      const pos_prec* qzold, const vel_prec* vlx,
                      const vel_prec* vly, const vel_prec* vlz);
void propagate_pos_lf_acc(time_prec dt, pos_prec* qx, pos_prec* qy,
                          pos_prec* qz, const pos_prec* qxold,
                          const pos_prec* qyold, const pos_prec* qzold,
                          const vel_prec* vlx, const vel_prec* vly,
                          const vel_prec* vlz);


void propagate_velocity_lp(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
                           const vel_prec* vxnew_lp, const vel_prec* vynew_lp,
                           const vel_prec* vznew_lp, const vel_prec* vxold_lp,
                           const vel_prec* vyold_lp, const vel_prec* vzold_lp,
                           const double scale, energy_prec& eksum_new,
                           energy_prec& eksum_old);
void propagate_velocity_lp_acc(
   vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp, const vel_prec* vxnew_lp,
   const vel_prec* vynew_lp, const vel_prec* vznew_lp, const vel_prec* vxold_lp,
   const vel_prec* vyold_lp, const vel_prec* vzold_lp, const double scale,
   energy_prec& eksum_new, energy_prec& eksum_old);


void propagate_velocity_lp2(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp,
                            vel_prec* vz_lp, const pos_prec* x_lp,
                            const pos_prec* y_lp, const pos_prec* z_lp,
                            const pos_prec* xold_lp, const pos_prec* yold_lp,
                            const pos_prec* zold_lp);
void propagate_velocity_lp2_acc(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp,
                                vel_prec* vz_lp, const pos_prec* x_lp,
                                const pos_prec* y_lp, const pos_prec* z_lp,
                                const pos_prec* xold_lp,
                                const pos_prec* yold_lp,
                                const pos_prec* zold_lp);


void propagate_velocity_lp3(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
                            const vel_prec* vxnew_lp, const vel_prec* vynew_lp,
                            const vel_prec* vznew_lp, const vel_prec* vxold_lp,
                            const vel_prec* vyold_lp, const vel_prec* vzold_lp,
                            energy_prec& eksum_new);
void propagate_velocity_lp3_acc(
   vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp, const vel_prec* vxnew_lp,
   const vel_prec* vynew_lp, const vel_prec* vznew_lp, const vel_prec* vxold_lp,
   const vel_prec* vyold_lp, const vel_prec* vzold_lp, energy_prec& eksum_new);
}
