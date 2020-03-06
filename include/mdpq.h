#pragma once
#include "md_prec.h"
#include "rc_man.h"
#include <istream>


/**
 * \defgroup md_pq  Atom Number, Momentum (p), and Coordinates (q)
 * \ingroup md
 */


TINKER_NAMESPACE_BEGIN
extern int rc_flag;


//====================================================================//


extern int n;
/**
 * \ingroup md_pq
 * \brief Number of atoms padded by #WARP_SIZE.
 * \see WARP_SIZE
 */
extern int padded_n;
/**
 * \ingroup md_pq
 * \brief Number of the trajectory frames.
 */
extern int trajn;


void n_data(rc_op);


//====================================================================//


/**
 * \ingroup md_pq
 * \brief Current coordinates used in energy evaluation and neighbor lists.
 */
extern real *x, *y, *z;
/**
 * \ingroup md_pq
 * \brief Entire trajectory frames.
 */
extern real *trajx, *trajy, *trajz;
/**
 * \ingroup md_pq
 * \brief Coordinates used in integrators.
 *
 * \note
 *    - New arrays will be allocated only if `sizeof(pos_prec) > sizeof(real)`,
 *    otherwise, they will be aliases of #x, #y, and #z.
 *    - Whenever #xpos, #ypos, #zpos get updated by integrators or barostats,
 *    #x, #y, #z must be updated immediately.
 *
 * \see pos_prec
 * \see real
 */
extern pos_prec *xpos, *ypos, *zpos;
static_assert(sizeof(pos_prec) >= sizeof(real),
              "Type pos_prec cannot be shorter than type real.");
/**
 * \ingroup md_pq
 * \brief Update #x, #y, #z by #xpos, #ypos, and #zpos.
 * If #xpos etc. are only aliaes, return directly.
 */
void copy_pos_to_xyz();
void copy_pos_to_xyz_acc();


/**
 * \ingroup md_pq
 * \brief Update #x, #y, #z via x += v * dt.
 * Currently #xpos, #ypos, and #zpos are integrated first, then call
 * #copy_pos_to_xyz() to update #x, #y, #z. In the end, neighbor lists may or
 * may not get updated.
 *
 * \param dt            Time-step for this update.
 * \param check_nblist  If `ture`, update the neighbor lists after updating the
 *                      coordinates.
 */
void propagate_xyz(time_prec dt, bool check_nblist);
void propagate_pos_acc(time_prec);


/**
 * \ingroup md_pq
 * \brief Finds the geometric center of each molecule and translate any stray
 * molecules back into the periodic box on GPU.
 * \note
 *    - Updating #x, #y, #z is the goal.
 *    - Checks whether PBC is in use inside the this function.
 *    - Will not perturb the neighbor lists so no need to update them.
 *    - Tinker uses center of mass.
 */
void bounds();
void bounds_pos_acc();
/**
 * \ingroup md_pq
 * \brief Call bounds() at least every x steps in MD.
 */
constexpr int BOUNDS_EVERY_X_STEPS = 500;

void read_frame_copyin_to_xyz(std::istream& input, int& done);


void xyz_data(rc_op);


//====================================================================//


/**
 * \ingroup md_pq
 * \brief Atomic mass.
 */
extern mass_prec* mass;
/**
 * \ingroup md_pq
 * \brief Inversed atomic mass.
 */
extern mass_prec* massinv;


/**
 * \ingroup md_pq
 * \brief Velocities.
 */
extern vel_prec *vx, *vy, *vz;


/**
 * \ingroup md_pq
 * \brief Update velocities via v += -g/m dt.
 */
void propagate_velocity(time_prec dt, const real* grx, const real* gry,
                        const real* grz);
/**
 * \ingroup md_pq
 * \brief Update velocities via v += -g/m dt.
 */
void propagate_velocity(time_prec dt, const fixed* grx, const fixed* gry,
                        const fixed* grz);
/**
 * \ingroup md_pq
 * \brief Update velocities via v += (-g/m dt -g2/m dt2).
 */
void propagate_velocity2(time_prec dt, const real* grx, const real* gry,
                         const real* grz, time_prec dt2, const real* grx2,
                         const real* gry2, const real* grz2);
/**
 * \ingroup md_pq
 * \brief Update velocities via v += (-g/m dt -g2/m dt2).
 */
void propagate_velocity2(time_prec dt, const fixed* grx, const fixed* gry,
                         const fixed* grz, time_prec dt2, const fixed* grx2,
                         const fixed* gry2, const fixed* grz2);
void propagate_velocity_acc(time_prec, const real*, const real*, const real*);
void propagate_velocity_acc(time_prec, const fixed*, const fixed*,
                            const fixed*);
void propagate_velocity2_acc(time_prec, const real*, const real*, const real*,
                             time_prec, const real*, const real*, const real*);
void propagate_velocity2_acc(time_prec, const fixed*, const fixed*,
                             const fixed*, time_prec, const fixed*,
                             const fixed*, const fixed*);


void mass_data(rc_op);
void vel_data(rc_op);
TINKER_NAMESPACE_END
