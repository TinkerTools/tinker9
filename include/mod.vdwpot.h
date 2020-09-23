#pragma once
#include "macro.h"
#include "mdprec.h"


namespace tinker {
/**
 * \ingroup vdw
 * \brief Constant flags for the VDW energy functions.
 */
enum class evdw_t
{
   decouple = 0,   ///< VDW lambda type: decouple.
   annihilate = 1, ///< VDW lambda type: annihilate.

   atom_type = 10,  ///< Indexing mode.
   atom_class = 11, ///< Indexing mode.

   arithmetic = 20, ///< Combining rule.
   geometric = 21,  ///< Combining rule.
   cubic_mean = 22, ///< Combining rule.
   hhg = 23,        ///< Combinding rule.

   lj,    ///< Lennard-Jones 12-6 potential.
   buck,  ///< Buckingham potential.
   mm3hb, ///< MM3 exp-6 potential.
   hal,   ///< Halgren buffered 14-7 potential.
   gauss, ///< Gaussian expansion VDW potential.
};
/**
 * \ingroup vdw
 * \brief Value of \f$ \gamma \f$ in buffered 14-7 vdw potential.
 */
TINKER_EXTERN real ghal;
/**
 * \ingroup vdw
 * \brief Value of \f$ \delta \f$ in buffered 14-7 vdw potential.
 */
TINKER_EXTERN real dhal;
TINKER_EXTERN real v2scale;
TINKER_EXTERN real v3scale;
TINKER_EXTERN real v4scale;
TINKER_EXTERN real v5scale;
TINKER_EXTERN evdw_t vdwtyp;
TINKER_EXTERN evdw_t vdwindex;
TINKER_EXTERN evdw_t radrule;
TINKER_EXTERN evdw_t epsrule;


//====================================================================//


/**
 * \ingroup vdw
 * \brief Long-range energy correction (lrc), used as `e += lrc/volume`.
 * \note Must be 0 if system is unbound.
 */
TINKER_EXTERN energy_prec elrc_vol;
/**
 * \ingroup vdw
 * \brief Long-range virial correction (lrc), used as `v(i,i) += lrc/volume`.
 * \note Must be 0 if system is unbound.
 */
TINKER_EXTERN virial_prec vlrc_vol;
}
