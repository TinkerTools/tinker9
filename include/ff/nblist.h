#pragma once
#include "precision.h"
#include "tool/genunit.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup nblist
enum class Nbl
{
   UNDEFINED = 0x00,   ///< Undefined.
   DOUBLE_LOOP = 0x01, ///< Double loop.
   VERLET = 0x02,      ///< Verlet neighbor list.
   SPATIAL = 0x04      ///< Spatial decomposition.
};
TINKER_ENABLE_ENUM_BITMASK(Nbl);

/// \ingroup nblist
/// \brief Verlet list: pairwise neighbor list indices and storage.
struct NBList
{
   int* nlst;     ///< number of sites in list for each atom
   int* lst;      ///< all of the sites in list
   int* update;   ///< update flag for each atom
   real* xold;    ///< old x coordinates
   real* yold;    ///< old y coordinates
   real* zold;    ///< old z coordinates
   const real* x; ///< current x coordinates
   const real* y; ///< current y coordinates
   const real* z; ///< current z coordinates
   int maxnlst;   ///< max number of neighbors for each atom
   real cutoff;   ///< list cutoff distance
   real buffer;   ///< width of the neighbor list buffer region

   ~NBList();
};
/// \ingroup nblist
using NBListUnit = GenericUnit<NBList, GenericUnitVersion::ENABLE_ON_DEVICE>;

/// \ingroup nblist
/// \brief
/// For Halgren Buffered 14-7 potential only, otherwise returns Nbl::UNDEFINED.
///
/// |                  | PBC                            | Unbound         |
/// |------------------|--------------------------------|-----------------|
/// | Inf. Cutoff      | N/A                            | double loop     |
/// | Cutoff + No List | double loop                    | double loop     |
/// | Cutoff + List    | verlet list or sptaial decomp. | verlet list (a) |
///
/// (a) We cannot use spatial decomposition because Tinker only set up a cubic
///     box for nonperiodic PME once. There is no guarantee that all of the atoms
///     will stay within the cubic box even if we moved the center of mass to
///     origin.
Nbl vlistVersion();
/// \ingroup nblist
/// \brief For partial charge models and for VDW models that do not
/// have a separate set of coordinates.
Nbl clistVersion();
/// \ingroup nblist
Nbl mlistVersion();
/// \ingroup nblist
Nbl ulistVersion();
/// \ingroup nblist
Nbl dsplistVersion();

/// \ingroup nblist
void nblistData(RcOp);

/// \ingroup nblist
void nblistRefresh();
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup nblist
/// \{
TINKER_EXTERN NBListUnit vlist_unit;
TINKER_EXTERN NBListUnit clist_unit;
TINKER_EXTERN NBListUnit mlist_unit;
TINKER_EXTERN NBListUnit ulist_unit;
TINKER_EXTERN NBListUnit dsplist_unit;
/// \}
}
