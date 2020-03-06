#pragma once
#include "box.h"
#include "darray.h"
#include "gen_unit.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
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
TINKER_EXTERN NBListUnit vlist_unit;
TINKER_EXTERN NBListUnit dlist_unit;
TINKER_EXTERN NBListUnit clist_unit;
TINKER_EXTERN NBListUnit mlist_unit;
TINKER_EXTERN NBListUnit ulist_unit;


/**
 * \brief
 * For Halgren Buffered 14-7 potential only, otherwise returns NBL_UNDEFINED.
 *
 */
nblist_t vlist_version();
nblist_t dlist_version();
nblist_t clist_version();
nblist_t mlist_version();
nblist_t ulist_version();


void nblist_data(rc_op op);


void nblist_build_acc(NBListUnit); // rc_init
void nblist_update_acc(NBListUnit);


void refresh_neighbors();
TINKER_NAMESPACE_END
