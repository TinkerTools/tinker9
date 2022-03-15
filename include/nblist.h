#pragma once
#include "box.h"
#include "tool/darray.h"
#include "tool/gen_unit.h"
#include "tool/rc_man.h"

namespace tinker {
enum class nblist_t
{
   UNDEFINED = 0x00,   ///< Undefined.
   DOUBLE_LOOP = 0x01, ///< Double loop.
   VERLET = 0x02,      ///< Verlet neighbor list.
   SPATIAL = 0x04      ///< Spatial decomposition.
};
TINKER_ENABLE_ENUM_BITMASK(nblist_t);
constexpr nblist_t NBL_UNDEFINED = nblist_t::UNDEFINED;
constexpr nblist_t NBL_DOUBLE_LOOP = nblist_t::DOUBLE_LOOP;
constexpr nblist_t NBL_VERLET = nblist_t::VERLET;
constexpr nblist_t NBL_SPATIAL = nblist_t::SPATIAL;

/// \brief
/// Verlet list: pairwise neighbor list indices and storage.
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
using NBListUnit = GenericUnit<NBList, GenericUnitVersion::EnableOnDevice>;

/**
 * \brief
 * For Halgren Buffered 14-7 potential only, otherwise returns NBL_UNDEFINED.
 *
 * |                  | PBC                            | Unbound         |
 * |------------------|--------------------------------|-----------------|
 * | Inf. Cutoff      | N/A                            | double loop     |
 * | Cutoff + No List | double loop                    | double loop     |
 * | Cutoff + List    | verlet list or sptaial decomp. | verlet list (a) |
 *
 * (a) We cannot use spatial decomposition because Tinker only set up a cubic
 *     box for nonperiodic PME once. There is no guarantee that all of the atoms
 *     will stay within the cubic box even if we moved the center of mass to
 *     origin.
 */
nblist_t vlist_version();
nblist_t dlist_version();
/**
 * \brief
 * For partial charge models and for VDW models that do not have a separate set
 * of coordinates.
 */
nblist_t clist_version();
nblist_t mlist_version();
nblist_t ulist_version();
nblist_t dsplist_version();

void nblist_data(rc_op op);

void nblist_build_acc(NBListUnit); // rc_init
void nblist_update_acc(NBListUnit);

void refresh_neighbors();
}
