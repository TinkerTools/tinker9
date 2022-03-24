#pragma once
#include "mod/elecamoeba.h"
#include "tool/rcman.h"

namespace tinker {
/**
 * \brief Multipole and AMOEBA polarization energy.
 * \note Will not be called in any of the following situations:
 *    - not using GPU;
 *    - not using CUDA as the primary GPU package;
 *    - not using periodic boundary condition;
 *    - not using both multipole and AMOEBA polarization terms;
 * \note Does not count number of interactions and aborts the program if called
 * erroneously (bug in the code).
 */
void emplar(int vers);
void emplar_cu(int);
}
