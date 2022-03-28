#include "ff/nblist.h"
#include "ff/amoeba/epolar.h"
#include "ff/elec.h"
#include "ff/hippo/edisp.h"
#include "ff/pchg/echarge.h"
#include "ff/pchg/echglj.h"
#include "ff/pchg/evdw.h"
#include "ff/potent.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "md/inc.h"
#include "ff/hippo/edisp.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/hippo/elechippo.h"
#include "ff/hippo/erepel.h"
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

Nbl vlistVersion()
{
   Nbl u;
   if (!usePotent(Potent::VDW)) {
      u = Nbl::UNDEFINED;
   } else if (vdwtyp != evdw_t::hal) {
      u = Nbl::UNDEFINED;
   } else if (!limits::use_vlist) {
      u = Nbl::DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = Nbl::VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = Nbl::SPATIAL;
      else
#endif
         u = Nbl::VERLET;
   }
   return u;
}

Nbl dlistVersion()
{
   return Nbl::UNDEFINED;
}

Nbl clistVersion()
{
   Nbl u;
   // First, forget about VDW, only check partial charge models.
   if (!usePotent(Potent::CHARGE) /* and !usePotent(Potent::SOLV) */) {
      u = Nbl::UNDEFINED;
   } else if (!limits::use_clist) {
      u = Nbl::DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = Nbl::VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = Nbl::SPATIAL;
      else
#endif
         u = Nbl::VERLET;
   }
   if (u != Nbl::UNDEFINED)
      return u;
   // Then, check VDW if no partial charge term is in use.
   if (!usePotent(Potent::VDW)) {
      u = Nbl::UNDEFINED;
   } else if (vdwtyp == evdw_t::hal) {
      u = Nbl::UNDEFINED;
   } else if (!limits::use_vlist) {
      u = Nbl::DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = Nbl::VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = Nbl::SPATIAL;
      else
#endif
         u = Nbl::VERLET;
   }
   return u;
}

Nbl mlistVersion()
{
   Nbl u;
   if (!usePotent(Potent::MPOLE) and !usePotent(Potent::POLAR) and !usePotent(Potent::CHGTRN) and
      !usePotent(Potent::REPULS) /* and !usePotent(Potent::SOLV) */) {
      u = Nbl::UNDEFINED;
   } else if (!limits::use_mlist) {
      u = Nbl::DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = Nbl::VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = Nbl::SPATIAL;
      else
#endif
         u = Nbl::VERLET;
   }
   return u;
}

Nbl ulistVersion()
{
   Nbl u;
   if (!usePotent(Potent::POLAR)) {
      u = Nbl::UNDEFINED;
   } else if (!limits::use_ulist) {
      u = Nbl::DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = Nbl::VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = Nbl::SPATIAL;
      else
#endif
         u = Nbl::VERLET;
   }
   return u;
}

Nbl dsplistVersion()
{
   Nbl u;
   if (!usePotent(Potent::DISP)) {
      u = Nbl::UNDEFINED;
   } else if (!limits::use_dlist) {
      u = Nbl::DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = Nbl::VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & Platform::CUDA)
         u = Nbl::SPATIAL;
      else
#endif
         u = Nbl::VERLET;
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
static void nblist_alloc(Nbl version, NBListUnit& nblu, int maxn, real cutoff, real buffer,
   const real* x, const real* y, const real* z)
{
   if (version & Nbl::DOUBLE_LOOP)
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
   SpatialUnit& unt, int n, real cut, real buf, const real* x, const real* y, const real* z,
   int nstype,                           //
   int ns1 = 0, int (*js1)[2] = nullptr, //
   int ns2 = 0, int (*js2)[2] = nullptr, //
   int ns3 = 0, int (*js3)[2] = nullptr, //
   int ns4 = 0, int (*js4)[2] = nullptr)
{
#if TINKER_CUDART
   spatialDataAlloc(unt, n, cut, buf, x, y, z, nstype, //
      ns1, js1, ns2, js2, ns3, js3, ns4, js4);
   alloc_thrust_cache = true;
#else
#endif
}

// rc_init
static void spatial_build(SpatialUnit unt)
{
#if TINKER_CUDART
   spatialDataInit_cu(unt);
#else
#endif
}

static void spatial_update(SpatialUnit unt)
{
#if TINKER_CUDART
   extern int check_spatial(
      int, real, int*, const real*, const real*, const real*, real*, real*, real*);
   auto& st = *unt;
   int answer =
      check_spatial(st.n, st.buffer, st.update, st.x, st.y, st.z, st.xold, st.yold, st.zold);
   if (answer) {
      spatialDataInit_cu(unt);
   } else {
      spatialDataUpdateSorted_cu(unt);
   }
#else
#endif
}

//====================================================================//

void nblist_build_acc(NBListUnit); // rc_init
void nblist_update_acc(NBListUnit);
void nblistData(RcOp op)
{
   if (op & rc_dealloc) {
      NBListUnit::clear();
      vlist_unit.close();
      clist_unit.close();
      mlist_unit.close();
      ulist_unit.close();
      dsplist_unit.close();

#if TINKER_CUDART
      SpatialUnit::clear();
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
   Nbl u = Nbl::UNDEFINED;
   double cut = 0;
   double buf = 0;

   // vlist
   u = vlistVersion();
   cut = switchOff(Switch::VDW);
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = vlist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, xred, yred, zred);
      }
      if (op & rc_init) {
         ehal_reduce_xyz();
         nblist_build_acc(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
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
   u = clistVersion();
   cut = -1;
   if (usePotent(Potent::CHARGE)) {
      cut = useEwald() ? switchOff(Switch::EWALD) : switchOff(Switch::CHARGE);
   }
   if (usePotent(Potent::VDW)) {
      double vdw_cut = switchOff(Switch::VDW);
      if (vdwtyp != evdw_t::hal)
         cut = std::max(cut, vdw_cut);
   }
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = clist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
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
   u = mlistVersion();
   cut = useEwald() ? switchOff(Switch::EWALD) : switchOff(Switch::MPOLE);
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = mlist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
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
   u = ulistVersion();
   cut = switchOff(Switch::USOLVE);
   buf = neigh::pbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = ulist_unit;
      if (op & rc_alloc) {
         const int maxnlst = 500;
         nblist_alloc(u, unt, maxnlst, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
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
   u = dsplistVersion();
   cut = useDEwald() ? switchOff(Switch::DEWALD) : switchOff(Switch::DISP);
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = dsplist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build_acc(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
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

void nblistRefresh()
{
   Nbl u = Nbl::UNDEFINED;

   // vlist
   u = vlistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = vlist_unit;
      ehal_reduce_xyz();
      nblist_update_acc(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = vspatial_v2_unit;
      ehal_reduce_xyz();
      spatial_update(un2);
   }

   // clist
   u = clistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = clist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = cspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }

   // mlist
   u = mlistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = mlist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = mspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }

   // ulist
   u = ulistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = ulist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = uspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatial_update(un2);
   }

   // dsplist
   u = dsplistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = dsplist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, g::q0);
      }
      nblist_update_acc(unt);
   }
   if (u & Nbl::SPATIAL) {
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
