#pragma once
#include "darray.h"
#include "rc_man.h"


namespace tinker {
TINKER_EXTERN int nmdpuexclude;
TINKER_EXTERN pointer<int, 2> mdpuexclude;
TINKER_EXTERN pointer<real, 4> mdpuexclude_scale;


/**
 * \note Directly return in any of the following situations:
 *    - not using GPU;
 *    - not using CUDA as the primary GPU package;
 *    - not using periodic boundary condition;
 *    - not using both multipole and AMOEBA polarization terms;
 * \note Shall not abort the program if called erroneously because of the early
 * return.
 */
void emplar_data(rc_op op);


/**
 * \brief Multipole and AMOEBA polarization energy.
 * \note Similar to emplar_data() but:
 *    - does not count number of interactions;
 *    - aborts the program if called erroneously (bug in the code).
 */
void emplar(int vers);
void emplar_cu(int);
}
