#pragma once
#include "ff/energybuffer.h"

// mpole
namespace tinker {
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
TINKER_EXTERN LocalFrame* zaxis;
TINKER_EXTERN real (*pole)[mpl_total];
TINKER_EXTERN real (*rpole)[mpl_total];
TINKER_EXTERN real* trqx;
TINKER_EXTERN real* trqy;
TINKER_EXTERN real* trqz;

TINKER_EXTERN virial_buffer vir_trq;

TINKER_EXTERN count_buffer nem;
TINKER_EXTERN energy_buffer em;
TINKER_EXTERN virial_buffer vir_em;
TINKER_EXTERN grad_prec* demx;
TINKER_EXTERN grad_prec* demy;
TINKER_EXTERN grad_prec* demz;
TINKER_EXTERN energy_prec energy_em;
TINKER_EXTERN virial_prec virial_em[9];
}
