#pragma once
#include "macro.h"


namespace tinker {
/**
 * \ingroup charge
 * \brief Magnitude of the partial charges (e-) of each **atom**.
 * \note Unlike Tinker, where only non-zero charges will be stored, this
 * array also includes zero charges.
 */
TINKER_EXTERN real* pchg;
}
