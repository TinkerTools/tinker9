#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
/**
 * \ingroup chgpen
 * \brief Charge Penetration damping type (GORDON1 or GORDON2)
 */
enum class chpdamp_t
{
   GORDON1,    ///< Lennard-Jones 12-6 potential.
   GORDON2,  ///< Buckingham potential.
};
TINKER_EXTERN real* pcore;
TINKER_EXTERN real* pval;
TINKER_EXTERN real* pval0;
TINKER_EXTERN real* palplha;


TINKER_EXTERN count_buffer ncp;


}
