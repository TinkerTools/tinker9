#pragma once
#include "energy_buffer.h"
#include "mod.energi.h"
#include "tool/rc_man.h"


namespace tinker {
extern real electric, dielec;
bool use_ewald();


//====================================================================//

/**
 * \ingroup charge
 * \brief Magnitude of the partial charges (e-) of each **atom**.
 * \note Unlike Tinker, where only non-zero charges will be stored, this
 * array also includes zero charges.
 */
extern real* pchg;
void pchg_data(rc_op);


//====================================================================//


enum
{
   mpl_pme_0 = 0,
   mpl_pme_x = 1,
   mpl_pme_y = 2,
   mpl_pme_z = 3,
   mpl_pme_xx = 4,
   mpl_pme_yy = 5,
   mpl_pme_zz = 6,
   mpl_pme_xy = 7,
   mpl_pme_xz = 8,
   mpl_pme_yz = 9,
   mpl_total = 10,
   mpl_pme_yx = mpl_pme_xy,
   mpl_pme_zx = mpl_pme_xz,
   mpl_pme_zy = mpl_pme_yz,


   pole_none = 0,
   pole_z_only = 1,
   pole_z_then_x = 2,
   pole_bisector = 3,
   pole_z_bisect = 4,
   pole_3_fold = 5
};


/**
 * \brief
 * Local axis type and x,y,z-axis defining atoms for each multipole site.
 */
struct LocalFrame
{
   int zaxis;  ///< Z-axis defining atom, starting from 0.
   int xaxis;  ///< X-axis defining atom, starting from 0.
   int yaxis;  ///< Y-axis defining atom, starting from ONE.
   int polaxe; ///< Local frame definition.
};
extern LocalFrame* zaxis;
extern real (*pole)[mpl_total];
extern real (*rpole)[mpl_total];
extern real *trqx, *trqy, *trqz;
extern real (*udir)[3];
extern real (*udirp)[3];
extern real (*uind)[3];
extern real (*uinp)[3];


void pole_data(rc_op);


//====================================================================//


void elec_data(rc_op op);


//====================================================================//


void mpole_init(int vers);
void chkpole();
void rotpole();
void torque(int vers);


void chkpole_acc();
void rotpole_acc();
void torque_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz);
// void torque_cu(int vers);
}
