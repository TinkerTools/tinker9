#pragma once

namespace tinker {
/// \ingroup vdw
/// \brief Constant flags for the VDW energy functions.
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
}
