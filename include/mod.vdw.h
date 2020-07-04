#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
/**
 * \ingroup vdw
 * \brief Type or class index into vdw parameters for each atom.
 * The indices have been sorted and start from 0.
 */
TINKER_EXTERN int* jvdw;
TINKER_EXTERN int* ired;
TINKER_EXTERN real* kred;
/**
 * \ingroup vdw
 * \brief Halgren buffered 14-7 reduced x, y, z coordinates for each atom.
 */
TINKER_EXTERN real* xred;
TINKER_EXTERN real* yred;
TINKER_EXTERN real* zred;
/**
 * \ingroup vdw
 * \brief Minimum energy distance (#radmin) or well depth parameter (#epsilon)
 * for each #jvdw pair. Element `[j1][j2]` is accessed by `[njvdw*j1 + j2]`.
 * \see njvdw
 */
TINKER_EXTERN real* radmin;
TINKER_EXTERN real* epsilon;
/**
 * \ingroup vdw
 * \brief VDW 1-4 parameters: minimum energy distance and well depth.
 * \see radmin epsilon
 */
TINKER_EXTERN real* radmin4;
TINKER_EXTERN real* epsilon4;


//====================================================================//


TINKER_EXTERN real* atom_rad;
TINKER_EXTERN real* atom_eps;


TINKER_EXTERN count_buffer nev;
TINKER_EXTERN energy_buffer ev;
TINKER_EXTERN virial_buffer vir_ev;
TINKER_EXTERN grad_prec* devx;
TINKER_EXTERN grad_prec* devy;
TINKER_EXTERN grad_prec* devz;
TINKER_EXTERN energy_prec energy_ev;
TINKER_EXTERN virial_prec virial_ev[9];


/**
 * \ingroup vdw
 * \brief Number of unique values in the #jvdw array.
 */
TINKER_EXTERN int njvdw;
/**
 * \ingroup vdw
 * \brief Halgren buffered 14-7 reduced vdw gradients for each atom.
 */
TINKER_EXTERN grad_prec* gxred;
TINKER_EXTERN grad_prec* gyred;
TINKER_EXTERN grad_prec* gzred;
TINKER_EXTERN bool vdw_exclude_bond;
TINKER_EXTERN int nvdw14;
TINKER_EXTERN int (*vdw14ik)[2];
TINKER_EXTERN int nvexclude;
TINKER_EXTERN int (*vexclude)[2];
TINKER_EXTERN real* vexclude_scale;
}
