#include "nblist.h"
#include "echarge.h"
#include "edisp.h"
#include "elec.h"
#include "epolar.h"
#include "evdw.h"
#include "glob.chglj.h"
#include "glob.nblist.h"
#include "glob.spatial.h"
#include "md.h"
#include "platform.h"
#include "potent.h"
#include "spatial2.h"
#include "switch.h"
#include "thrust_cache.h"
#include "tool/darray.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/limits.hh>
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
      if (pltfm_config & CU_PLTFM)
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
   if (!use_potent(charge_term) /* && !use_potent(solv_term) */) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_clist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & CU_PLTFM)
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
      if (pltfm_config & CU_PLTFM)
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
   if (!use_potent(mpole_term) && !use_potent(polar_term) &&
       !use_potent(chgtrn_term) &&
       !use_potent(repuls_term) /* && !use_potent(solv_term) */) {
      u = NBL_UNDEFINED;
   } else if (!limits::use_mlist) {
      u = NBL_DOUBLE_LOOP;
   } else if (!bound::use_bounds) {
      u = NBL_VERLET;
   } else {
#if TINKER_CUDART
      if (pltfm_config & CU_PLTFM)
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
      if (pltfm_config & CU_PLTFM)
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
      if (pltfm_config & CU_PLTFM)
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
static void nblist_alloc(nblist_t version, NBListUnit& nblu, int maxn,
                         real cutoff, real buffer, const real* x, const real* y,
                         const real* z)
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

   nblu.update_deviceptr(st, WAIT_NEW_Q);
}


#if TINKER_CUDART
static bool alloc_thrust_cache;
// rc_alloc
static void spatial_alloc(SpatialUnit& unt, int n, real cut, real buf,
                          const real* x, const real* y, const real* z)
{
   spatial_data_alloc(unt, n, cut, buf, x, y, z);
   alloc_thrust_cache = true;
}


// rc_init
template <class SPT>
static void spatial_build(SPT unt)
{
   spatial_data_init_cu(unt);
}


template <class SPT>
static void spatial_update(SPT unt)
{
   extern int check_spatial(int, real, int*, const real*, const real*,
                            const real*, real*, real*, real*);
   auto& st = *unt;
   int answer = check_spatial(st.n, st.buffer, st.update, st.x, st.y, st.z,
                              st.xold, st.yold, st.zold);
   if (answer) {
      spatial_data_init_cu(unt);
   } else {
      spatial_data_update_sorted(unt);
   }
}
#else
static void spatial_alloc(SpatialUnit&, int, real, real, const real*,
                          const real*, const real*)
{}
template <class SPT>
static void spatial_build(SPT)
{}
template <class SPT>
static void spatial_update(SPT)
{}
#endif
static void spatial_alloc( //
   Spatial2Unit& unt, int n, real cut, real buf, const real* x, const real* y,
   const real* z, int nstype,            //
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


//====================================================================//


void nblist_data(rc_op op)
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
      thrust_cache_dealloc();
      mspatial_unit.close();


      Spatial2Unit::clear();
      cspatial_v2_unit.close();
      vspatial_v2_unit.close();
      uspatial_v2_unit.close();
      dspspatial_v2_unit.close();
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
   cut = switch_off(switch_vdw);
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
         spatial_alloc(un2, n, cut, buf, xred, yred, zred, 1, nvexclude,
                       vexclude);
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
      cut = use_ewald() ? switch_off(switch_ewald) : switch_off(switch_charge);
   }
   if (use_potent(vdw_term)) {
      double vdw_cut = switch_off(switch_vdw);
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
         spatial_alloc(un2, n, cut, buf, x, y, z, 2, //
                       ncexclude, cexclude, nvexclude, vexclude);
      }
      if (op & rc_init) {
         spatial_build(un2);
      }
   }


   // mlist
   u = mlist_version();
   cut = use_ewald() ? switch_off(switch_ewald) : switch_off(switch_mpole);
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
      auto& unt = mspatial_unit;
      if (op & rc_alloc) {
         spatial_alloc(unt, n, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         spatial_build(unt);
      }
   }


   // ulist
   u = ulist_version();
   cut = switch_off(switch_usolve);
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
         spatial_alloc(un2, n, cut, buf, x, y, z, 1, nuexclude, uexclude);
      }
      if (op & rc_init) {
         spatial_build(un2);
      }
   }


   // dsplist
   u = dsplist_version();
   cut = use_dewald() ? switch_off(switch_dewald) : switch_off(switch_disp);
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
      thrust_cache_alloc();
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
         unt.update_deviceptr(*unt, PROCEED_NEW_Q);
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
         unt.update_deviceptr(*unt, PROCEED_NEW_Q);
      }
      nblist_update_acc(unt);
   }
   if (u & NBL_SPATIAL) {
      auto& unt = mspatial_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, PROCEED_NEW_Q);
      }
      spatial_update(unt);
   }


   // ulist
   u = ulist_version();
   if (u & (NBL_DOUBLE_LOOP | NBL_VERLET)) {
      auto& unt = ulist_unit;
      if (rc_flag & calc::traj) {
         unt->x = x;
         unt->y = y;
         unt->z = z;
         unt.update_deviceptr(*unt, PROCEED_NEW_Q);
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
         unt.update_deviceptr(*unt, PROCEED_NEW_Q);
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
