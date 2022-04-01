#include "ff/nblist.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/atom.h"
#include "ff/elec.h"
#include "ff/hippo/edisp.h"
#include "ff/hippo/elechippo.h"
#include "ff/hippo/erepel.h"
#include "ff/pchg/echarge.h"
#include "ff/pchg/echglj.h"
#include "ff/pchg/evdw.h"
#include "ff/potent.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "tool/thrustcache.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/neigh.hh>

namespace tinker {
NBList::~NBList()
{
   darray::deallocate(nlst, lst, update, xold, yold, zold);
}
}

namespace tinker {
Nbl vlistVersion()
{
   Nbl u;
   if (not usePotent(Potent::VDW)) {
      u = Nbl::UNDEFINED;
   } else if (vdwtyp != Vdw::HAL) {
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

Nbl clistVersion()
{
   Nbl u;
   // First, forget about VDW, only check partial charge models.
   if (not usePotent(Potent::CHARGE) /* and not usePotent(Potent::SOLV) */) {
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
   if (not usePotent(Potent::VDW)) {
      u = Nbl::UNDEFINED;
   } else if (vdwtyp == Vdw::HAL) {
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
   if (not usePotent(Potent::MPOLE) and not usePotent(Potent::POLAR) and
      not usePotent(Potent::CHGTRN) and
      not usePotent(Potent::REPULS) /* and not usePotent(Potent::SOLV) */) {
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
   if (not usePotent(Potent::POLAR)) {
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
   if (not usePotent(Potent::DISP)) {
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
}

namespace tinker {
#if TINKER_CUDART
static bool alloc_thrust_cache;
#endif

// In the gas phase calculation where neighbor list is not used, we should
// always first check the value of `maxn`. If `maxn` is equal to 1, it means
// the value of cutoff can even be `INF`.
// \see cutoffs.f
static int nblistMaxlst(int maxn, double cutoff, double buffer)
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

// RcOp::ALLOC
static void nblistAlloc(Nbl version, NBListUnit& nblu, int maxn, real cutoff, real buffer,
   const real* x, const real* y, const real* z)
{
   if (version & Nbl::DOUBLE_LOOP)
      maxn = 1;

   nblu = NBListUnit::open();
   auto& st = *nblu;

   darray::allocate(n, &st.nlst);

   int maxlst = nblistMaxlst(maxn, cutoff, buffer);
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

   nblu.deviceptrUpdate(st, g::q0);
   waitFor(g::q0);
}

// RcOp::ALLOC
static void spatialAlloc( //
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

// RcOp::INIT
extern void nblistBuild_acc(NBListUnit);
static void nblistBuild(NBListUnit u)
{
   nblistBuild_acc(u);
}

// RcOp::INIT
static void spatialBuild(SpatialUnit unt)
{
#if TINKER_CUDART
   spatialDataInit(unt);
#else
#endif
}

void nblistData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
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

   if (op & RcOp::ALLOC) {
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
      if (op & RcOp::ALLOC) {
         nblistAlloc(u, unt, 2500, cut, buf, xred, yred, zred);
      }
      if (op & RcOp::INIT) {
         ehalReduceXyz();
         nblistBuild(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = vspatial_v2_unit;
      if (op & RcOp::ALLOC) {
         spatialAlloc(un2, n, cut, buf, xred, yred, zred, 1, nvexclude, vexclude);
      }
      if (op & RcOp::INIT) {
         ehalReduceXyz();
         spatialBuild(un2);
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
      if (vdwtyp != Vdw::HAL)
         cut = std::max(cut, vdw_cut);
   }
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = clist_unit;
      if (op & RcOp::ALLOC) {
         nblistAlloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & RcOp::INIT) {
         nblistBuild(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = cspatial_v2_unit;
      if (op & RcOp::ALLOC) {
         spatialAlloc(un2, n, cut, buf, x, y, z, 3, //
            ncexclude, cexclude, nvexclude, vexclude, ncvexclude, cvexclude);
      }
      if (op & RcOp::INIT) {
         spatialBuild(un2);
      }
   }

   // mlist
   u = mlistVersion();
   cut = useEwald() ? switchOff(Switch::EWALD) : switchOff(Switch::MPOLE);
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = mlist_unit;
      if (op & RcOp::ALLOC) {
         nblistAlloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & RcOp::INIT) {
         nblistBuild(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = mspatial_v2_unit;
      if (op & RcOp::ALLOC) {
         if (mplpot::use_chgpen) {
            spatialAlloc(
               un2, n, cut, buf, x, y, z, 2, nmdwexclude, mdwexclude, nrepexclude, repexclude);
         } else {
            spatialAlloc(un2, n, cut, buf, x, y, z, 4, nmdpuexclude, mdpuexclude, nmexclude,
               mexclude, ndpexclude, dpexclude, nuexclude, uexclude);
         }
      }
      if (op & RcOp::INIT) {
         spatialBuild(un2);
      }
   }

   // ulist
   u = ulistVersion();
   cut = switchOff(Switch::USOLVE);
   buf = neigh::pbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = ulist_unit;
      if (op & RcOp::ALLOC) {
         const int maxnlst = 500;
         nblistAlloc(u, unt, maxnlst, cut, buf, x, y, z);
      }
      if (op & RcOp::INIT) {
         nblistBuild(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = uspatial_v2_unit;
      if (op & RcOp::ALLOC) {
         if (mplpot::use_chgpen) {
            spatialAlloc(un2, n, cut, buf, x, y, z, 1, nwexclude, wexclude);
         } else {
            spatialAlloc(un2, n, cut, buf, x, y, z, 1, nuexclude, uexclude);
         }
      }
      if (op & RcOp::INIT) {
         spatialBuild(un2);
      }
   }

   // dsplist
   u = dsplistVersion();
   cut = useDEwald() ? switchOff(Switch::DEWALD) : switchOff(Switch::DISP);
   buf = neigh::lbuffer;
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = dsplist_unit;
      if (op & RcOp::ALLOC) {
         nblistAlloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & RcOp::INIT) {
         nblistBuild(unt);
      }
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = dspspatial_v2_unit;
      if (op & RcOp::ALLOC) {
         spatialAlloc(un2, n, cut, buf, x, y, z, 1, ndspexclude, dspexclude);
      }
      if (op & RcOp::INIT) {
         spatialBuild(un2);
      }
   }

#if TINKER_CUDART
   if (alloc_thrust_cache)
      ThrustCache::allocate();
#endif
}
}

namespace tinker {
extern void nblistUpdate_acc(NBListUnit);
static void nblistUpdate(NBListUnit u)
{
   nblistUpdate_acc(u);
}

extern int spatialCheck_acc(
   int, real, int*, const real*, const real*, const real*, real*, real*, real*);
static void spatialUpdate(SpatialUnit unt)
{
#if TINKER_CUDART
   auto& st = *unt;
   int answer =
      spatialCheck_acc(st.n, st.buffer, st.update, st.x, st.y, st.z, st.xold, st.yold, st.zold);
   if (answer) {
      spatialDataInit(unt);
   } else {
      spatialDataUpdateSorted(unt);
   }
#endif
}

void nblistRefresh()
{
   Nbl u = Nbl::UNDEFINED;

   // vlist
   u = vlistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = vlist_unit;
      ehalReduceXyz();
      nblistUpdate(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = vspatial_v2_unit;
      ehalReduceXyz();
      spatialUpdate(un2);
   }

   // clist
   u = clistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = clist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.deviceptrUpdate(*unt, g::q0);
      }
      nblistUpdate(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = cspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatialUpdate(un2);
   }

   // mlist
   u = mlistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = mlist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.deviceptrUpdate(*unt, g::q0);
      }
      nblistUpdate(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = mspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatialUpdate(un2);
   }

   // ulist
   u = ulistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = ulist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.deviceptrUpdate(*unt, g::q0);
      }
      nblistUpdate(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = uspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatialUpdate(un2);
   }

   // dsplist
   u = dsplistVersion();
   if (u & (Nbl::DOUBLE_LOOP | Nbl::VERLET)) {
      auto& unt = dsplist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.deviceptrUpdate(*unt, g::q0);
      }
      nblistUpdate(unt);
   }
   if (u & Nbl::SPATIAL) {
      auto& un2 = dspspatial_v2_unit;
      if (rc_flag & calc::traj) {
         un2->x = x;
         un2->y = y;
         un2->z = z;
      }
      spatialUpdate(un2);
   }
}
}
