#include "ff/nblist.h"
#include "ff/amoeba/epolar.h"
#include "ff/elec.h"
#include "ff/hippo/edisp.h"
#include "ff/pchg/echarge.h"
#include "ff/pchg/evdw.h"
#include "ff/potent.h"
#include "ff/spatial2.h"
#include "ff/switch.h"
#include "md/inc.h"
#include "mod/disp.h"
#include "mod/elecamoeba.h"
#include "mod/elechippo.h"
#include "mod/nblist.h"
#include "mod/repel.h"
#include "platform.h"
#include "tool/darray.h"
#include "tool/thrustcache.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/neigh.hh>
#include <tinker/detail/vdwpot.hh>

namespace tinker {
NBList::~NBList()
{
   darray::deallocate(nlst, lst, update, xold, yold, zold);
}

//====================================================================//

nblist_t vlist_version()
{
   nblist_t u;
   if (!use_potent(vdw_term)) {
      u = NBL_UNDEFINED;
   } else if (vdwtyp != evdw_t::hal) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_vlist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = NBL_SPATIAL;
      else
#endif
         u = NBL_VERLET;
   }
   return u;
}

nblist_t dlist_version()
{
   return NBL_UNDEFINED;
}

nblist_t clist_version()
{
   nblist_t u;
   // First, forget about VDW, only check partial charge models.
   if (!use_potent(charge_term) /* and !use_potent(solv_term) */) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_clist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = NBL_SPATIAL;
      else
#endif
         u = NBL_VERLET;
   }
   if (u != NBL_UNDEFINED)
      return u;
   // Then, check VDW if no partial charge term is in use.
   if (!use_potent(vdw_term)) {
      u = NBL_UNDEFINED;
   } else if (vdwtyp == evdw_t::hal) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_vlist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = NBL_SPATIAL;
      else
#endif
         u = NBL_VERLET;
   }
   return u;
}

nblist_t mlist_version()
{
   nblist_t u;
   if (!use_potent(mpole_term) and !use_potent(polar_term) and !use_potent(chgtrn_term) and
      !use_potent(repuls_term) /* and !use_potent(solv_term) */) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_mlist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = NBL_SPATIAL;
      else
#endif
         u = NBL_VERLET;
   }
   return u;
}

nblist_t ulist_version()
{
   nblist_t u;
   if (!use_potent(polar_term)) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_ulist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = NBL_SPATIAL;
      else
#endif
         u = NBL_VERLET;
   }
   return u;
}

nblist_t dsplist_version()
{
   nblist_t u;
   if (!use_potent(disp_term)) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_dlist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = NBL_SPATIAL;
      else
#endif
         u = NBL_VERLET;
   }
   return u;
}

//====================================================================//

/**
 * In the gas phase calculation where neighbor list is not used, we should
 * always first check the value of `maxn`. If `maxn` is equal to 1, it means
 * the value of cutoff can even be `INF`.
 * \see cutoffs.f
 */
static int nblist_maxlst(int maxn, double cutoff, double buffer)
{
   if (maxn > 1) {
      double buf = (cutoff + buffer);
      double buf3 = buf * buf * buf + 100; // empirical formula
      int limit;
      // assuming buf3 does not overflow
      // compare buf3 to 0x10000000 while max of int == 0x7FFFFFFF
      if (buf3 > 0x10000000)
         limit = maxn;
      else
         limit = buf3;
      int ans = std::min(limit, maxn);
      if (ans > 1) {
         const int magic = 32;
         ans = (ans + magic - 1) / magic;
         ans *= magic;
      }
      return ans;
   } else {
      return 1;
   }
}

// rc_alloc
static void nblist_alloc(nblist_t version, NBListUnit& nblu, int maxn, real cutoff, real buffer,
   const real* x, const real* y, const real* z)
{
   if (version & NBL_DOUBLE_LOOP)
      maxn = 1;

   nblu = NBListUnit::open();
   auto& st = *nblu;

   darray::allocate(n, &st.nlst);

   int maxlst = nblist_maxlst(maxn, cutoff, buffer);
   darray::allocate(maxlst * n, &st.lst);

   if (maxlst == 1) {
      st.update = nullptr;
      st.xold = nullptr;
      st.yold = nullptr;
      st.zold = nullptr;
   } else {
      darray::allocate(n, &st.update, &st.xold, &st.yold, &st.zold);
   }

   st.x = x;
   st.y = y;
   st.z = z;

   st.maxnlst = maxlst;
   st.cutoff = cutoff;
   st.buffer = buffer;

   nblu.update_deviceptr(st, g::q0);
   wait_for(g::q0);
}

#if TINKER_CUDART
static bool alloc_thrust_cache;
#else
#endif

// rc_alloc
static void spatial_alloc( //
   Spatial2Unit& unt, int n, real cut, real buf, const real* x, const real* y, const real* z,
   int nstype,                           //
   int ns1 = 0, int (*js1)[2] = nullptr, //
   int ns2 = 0, int (*js2)[2] = nullptr, //
   int ns3 = 0, int (*js3)[2] = nullptr, //
   int ns4 = 0, int (*js4)[2] = nullptr)
{
#if TINKER_CUDART
   spatial2_data_alloc(unt, n, cut, buf, x, y, z, nstype, //
      ns1, js1, ns2, js2, ns3, js3, ns4, js4);
   alloc_thrust_cache = true;
#else
#endif
}

// rc_init
static void spatial_build(Spatial2Unit unt)
{
#if TINKER_CUDART
   spatial_data_init_cu(unt);
#else
#endif
}

static void spatial_update(Spatial2Unit unt)
{
#if TINKER_CUDART
   extern int check_spatial(
      int, real, int*, const real*, const real*, const real*, real*, real*, real*);
   auto& st = *unt;
   int answer =
      check_spatial(st.n, st.buffer, st.update, st.x, st.y, st.z, st.xold, st.yold, st.zold);
   if (answer) {
      spatial_data_init_cu(unt);
   } else {
      spatial_data_update_sorted(unt);
   }
#else
#endif
}

//====================================================================//

void nblist_data(RcOp op)
{
   if (op & rc_dealloc) {
      NBListUnit::clear();
      vlist_unit.close();
      clist_unit.close();
      mlist_unit.close();
      ulist_unit.close();
      dsplist_unit.close();

#if TINKER_CUDART
      Spatial2Unit::clear();
      cspatial_v2_unit.close();
      vspatial_v2_unit.close();
      uspatial_v2_unit.close();
      mspatial_v2_unit.close();
      dspspatial_v2_unit.close();

      ThrustCache::deallocate();
#endif
   }

   if (op & rc_alloc) {
      assert(NBListUnit::size() == 0);
      assert(SpatialUnit::size() == 0);
   }

#if TINKER_CUDART
   alloc_thrust_cache = false;
#endif
   nblist_t u = NBL_UNDEFINED;
   double cut = 0;
   double buf = 0;

   // vlist
   u = vlist_version();
   cut = switchOff(SWITCH_VDW);
   buf = neigh::lbuffer;
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = vlist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, xred, yred, zred);
      }
      if (op & rc_init) {
         ehal_reduce_xyz();
         nblist_build_acc(unt);
      }
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = vspatial_v2_unit;
      if (op & rc_alloc) {
         spatial_alloc(un2, n, cut, buf, xred, yred, zred, 1, nvexclude, vexclude);
      }
      if (op & rc_init) {
         ehal_reduce_xyz();
         spatial_build(un2);
      }
   }

   // clist
   u = clist_version();
   cut = -1;
   if (use_potent(charge_term)) {
      cut = useEwald() ? switchOff(SWITCH_EWALD) : switchOff(SWITCH_CHARGE);
   }
   if (use_potent(vdw_term)) {
      double vdw_cut = switchOff(SWITCH_VDW);
      if (vdwtyp != evdw_t::hal)
         cut = std::max(cut, vdw_cut);
   }
   buf = neigh::lbuffer;
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = clist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = cspatial_v2_unit;
      if (op & rc_alloc) {
         spatial_alloc(un2, n, cut, buf, x, y, z, 3, //
            ncexclude, cexclude, nvexclude, vexclude, ncvexclude, cvexclude);
      }
      if (op & rc_init) {
         spatial_build(un2);
      }
   }

   // mlist
   u = mlist_version();
   cut = useEwald() ? switchOff(SWITCH_EWALD) : switchOff(SWITCH_MPOLE);
   buf = neigh::lbuffer;
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = mlist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = mspatial_v2_unit;
      if (op & rc_alloc) {
         if (mplpot::use_chgpen) {
            spatial_alloc(
               un2, n, cut, buf, x, y, z, 2, nmdwexclude, mdwexclude, nrepexclude, repexclude);
         } else {
            spatial_alloc(un2, n, cut, buf, x, y, z, 4, nmdpuexclude, mdpuexclude, nmexclude,
               mexclude, ndpexclude, dpexclude, nuexclude, uexclude);
         }
      }
      if (op & rc_init) {
         spatial_build(un2);
      }
   }

   // ulist
   u = ulist_version();
   cut = switchOff(SWITCH_USOLVE);
   buf = neigh::pbuffer;
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = ulist_unit;
      if (op & rc_alloc) {
         const int maxnlst = 500;
         nblist_alloc(u, unt, maxnlst, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = uspatial_v2_unit;
      if (op & rc_alloc) {
         if (mplpot::use_chgpen) {
            spatial_alloc(un2, n, cut, buf, x, y, z, 1, nwexclude, wexclude);
         } else {
            spatial_alloc(un2, n, cut, buf, x, y, z, 1, nuexclude, uexclude);
         }
      }
      if (op & rc_init) {
         spatial_build(un2);
      }
   }

   // dsplist
   u = dsplist_version();
   cut = useDEwald() ? switchOff(SWITCH_DEWALD) : switchOff(SWITCH_DISP);
   buf = neigh::lbuffer;
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = dsplist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = dspspatial_v2_unit;
      if (op & rc_alloc) {
         spatial_alloc(un2, n, cut, buf, x, y, z, 1, ndspexclude, dspexclude);
      }
      if (op & rc_init) {
         spatial_build(un2);
      }
   }

#if TINKER_CUDART
   if (alloc_thrust_cache)
      ThrustCache::allocate();
#endif
}

void refresh_neighbors()
{
   nblist_t u = NBL_UNDEFINED;

   // vlist
   u = vlist_version();
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = vlist_unit;
      ehal_reduce_xyz();
      nblist_update_acc(unt);
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = vspatial_v2_unit;
      ehal_reduce_xyz();
      spatial_update(un2);
   }

   // clist
   u = clist_version();
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = clist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = cspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }

   // mlist
   u = mlist_version();
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = mlist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = mspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }

   // ulist
   u = ulist_version();
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = ulist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = uspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }

   // dsplist
   u = dsplist_version();
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = dsplist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & NBL_SPATIAL) {
      auto& un2 = dspspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }
}
}
