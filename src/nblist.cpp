#include "nblist.h"
#include "dev_array.h"
#include "e_polar.h"
#include "e_vdw.h"
#include "elec.h"
#include "md.h"
#include "potent.h"
#include "spatial.h"
#include "switch.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/neigh.hh>


TINKER_NAMESPACE_BEGIN
int always_use_nblist;


NBList::~NBList()
{
   device_array::deallocate(nlst, lst, update, xold, yold, zold);
}


//====================================================================//


int vlist_version()
{
   if (!use_potent(vdw_term))
      return 0;
   if (!limits::use_vlist)
      return NBList::double_loop;
   if (!bound::use_bounds)
      return NBList::nblist;
#if TINKER_CUDART
   if (always_use_nblist)
      return NBList::nblist;
   else
      return NBList::spatial;
#else
   return NBList::nblist;
#endif
}


int dlist_version()
{
   if (!use_potent(disp_term))
      return 0;
   if (!limits::use_dlist)
      return NBList::double_loop;
   if (!bound::use_bounds)
      return NBList::nblist;
#if TINKER_CUDART
   return NBList::spatial;
#else
   return NBList::nblist;
#endif
}


int clist_version()
{
   if (!use_potent(charge_term) /* && !use_potent(solv_term) */)
      return 0;
   if (!limits::use_clist)
      return NBList::double_loop;
   if (!bound::use_bounds)
      return NBList::nblist;
#if TINKER_CUDART
   return NBList::spatial;
#else
   return NBList::nblist;
#endif
}


int mlist_version()
{
   if (!use_potent(mpole_term) && !use_potent(polar_term)
       /* && !use_potent(chgtrn_term) && !use_potent(solv_term) */)
      return 0;
   if (!limits::use_mlist)
      return NBList::double_loop;
   if (!bound::use_bounds)
      return NBList::nblist;
#if TINKER_CUDART
   return NBList::spatial;
#else
   return NBList::nblist;
#endif
}


int ulist_version()
{
   if (!use_potent(polar_term))
      return 0;
   if (!limits::use_ulist)
      return NBList::double_loop;
   if (!bound::use_bounds)
      return NBList::nblist;
#if TINKER_CUDART
   return NBList::spatial;
#else
   return NBList::nblist;
#endif
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
static void nblist_alloc(int version, NBListUnit& nblu, int maxn, real cutoff,
                         real buffer, const real* x, const real* y,
                         const real* z)
{
   if (version & NBList::double_loop)
      maxn = 1;


   nblu = NBListUnit::open();
   auto& st = *nblu;

   device_array::allocate(n, &st.nlst);

   int maxlst = nblist_maxlst(maxn, cutoff, buffer);
   device_array::allocate(maxlst * n, &st.lst);

   if (maxlst == 1) {
      st.update = nullptr;
      st.xold = nullptr;
      st.yold = nullptr;
      st.zold = nullptr;
   } else {
      device_array::allocate(n, &st.update, &st.xold, &st.yold, &st.zold);
   }

   st.x = x;
   st.y = y;
   st.z = z;

   st.maxnlst = maxlst;
   st.cutoff = cutoff;
   st.buffer = buffer;

   nblu.update_deviceptr(st);
}


// rc_init
extern void nblist_build(NBListUnit);


// rc_evolve
extern void nblist_update(NBListUnit);


static bool alloc_thrust_cache;
extern void spatial_data_init_cu(SpatialUnit);
#if TINKER_CUDART
// rc_alloc
static void spatial_alloc(SpatialUnit& unt, int n, real cut, real buf,
                          const real* x, const real* y, const real* z)
{
   spatial_data_alloc(unt, n, cut, buf, x, y, z);
   alloc_thrust_cache = true;
}


// rc_init
static void spatial_build(SpatialUnit unt)
{
   spatial_data_init_cu(unt);
}


// rc_evolve
static void spatial_update(SpatialUnit unt)
{
   extern int check_spatial(int, real, const Box*, int*, const real*,
                            const real*, const real*, real*, real*, real*);
   auto& st = *unt;
   st.rebuild = check_spatial(st.n, st.buffer, box, st.update, st.x, st.y, st.z,
                              st.xold, st.yold, st.zold);
   if (st.rebuild)
      spatial_data_init_cu(unt);
}
#else
static void spatial_alloc(SpatialUnit& unt, int n, real cutoff, real buffer,
                          const real* x, const real* y, const real* z)
{}
static void spatial_build(SpatialUnit) {}
static void spatial_update(SpatialUnit) {}
#endif


//====================================================================//


void nblist_data(rc_op op)
{
   if (op & rc_dealloc) {
      NBListUnit::clear();
      vlist_unit.close();
      dlist_unit.close();
      clist_unit.close();
      mlist_unit.close();
      ulist_unit.close();

      always_use_nblist = 0;

#if TINKER_CUDART
      extern void thrust_cache_dealloc();
      SpatialUnit::clear();
      thrust_cache_dealloc();
      vspatial_unit.close();
      mspatial_unit.close();
#endif
   }


   if (op & rc_alloc) {
      assert(NBListUnit::size() == 0);
      assert(SpatialUnit::size() == 0);
   }


   alloc_thrust_cache = false;
   int u = 0;
   double cut = 0;
   double buf = 0;


   // vlist
   u = vlist_version();
   cut = switch_off(switch_vdw);
   buf = neigh::lbuffer;
   if (u & (NBList::double_loop | NBList::nblist)) {
      auto& unt = vlist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, xred, yred, zred);
      }
      if (op & rc_init) {
         evdw_reduce_xyz();
         nblist_build(unt);
      }
      if (op & rc_evolve) {
         evdw_reduce_xyz();
         nblist_update(unt);
      }
   }
   if (u & NBList::spatial) {
      auto& unt = vspatial_unit;
      if (op & rc_alloc) {
         spatial_alloc(unt, n, cut, buf, xred, yred, zred);
      }
      if (op & rc_init) {
         evdw_reduce_xyz();
         spatial_build(unt);
      }
      if (op & rc_evolve) {
         evdw_reduce_xyz();
         spatial_update(unt);
      }
   }


   // dlist


   // clist


   // mlist
   u = mlist_version();
   cut = use_ewald() ? switch_off(switch_ewald) : switch_off(switch_mpole);
   buf = neigh::lbuffer;
   if (u & (NBList::double_loop | NBList::nblist)) {
      auto& unt = mlist_unit;
      if (op & rc_alloc) {
         nblist_alloc(u, unt, 2500, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         nblist_build(unt);
      }
      if (op & rc_evolve) {
         if (rc_flag & calc::traj) {
            unt->x = x;
            unt->y = y;
            unt->z = z;
            unt.update_deviceptr(*unt);
         }
         nblist_update(unt);
      }
   }
   if (u & NBList::spatial) {
      auto& unt = mspatial_unit;
      if (op & rc_alloc) {
         spatial_alloc(unt, n, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         spatial_build(unt);
      }
      if (op & rc_evolve) {
         if (rc_flag & calc::traj) {
            unt->x = x;
            unt->y = y;
            unt->z = z;
            unt.update_deviceptr(*unt);
         }
         spatial_update(unt);
      }
   }


   // ulist
   u = ulist_version();
   cut = switch_off(switch_usolve);
   buf = neigh::pbuffer;
   if (u & (NBList::double_loop | NBList::nblist)) {
      auto& unt = ulist_unit;
      if (op & rc_alloc) {
         const int maxnlst = 500;
         nblist_alloc(u, unt, maxnlst, cut, buf, x, y, z);
         device_array::allocate(n, &mindex);
         int minv_size = nblist_maxlst(maxnlst, cut, buf);
         device_array::allocate(3 * minv_size * n, &minv);
         device_array::allocate(6 * nuexclude_, &minv_exclude_);
      }
      if (op & rc_init) {
         nblist_build(unt);
      }
      if (op & rc_evolve) {
         if (rc_flag & calc::traj) {
            unt->x = x;
            unt->y = y;
            unt->z = z;
            unt.update_deviceptr(*unt);
         }
         nblist_update(unt);
      }
      if (op & rc_dealloc) {
         device_array::deallocate(mindex, minv, minv_exclude_);
      }
   }
   if (u & NBList::spatial) {
      auto& unt = uspatial_unit;
      if (op & rc_alloc) {
         spatial_alloc(unt, n, cut, buf, x, y, z);
      }
      if (op & rc_init) {
         spatial_build(unt);
      }
      if (op & rc_evolve) {
         if (rc_flag & calc::traj) {
            unt->x = x;
            unt->y = y;
            unt->z = z;
            unt.update_deviceptr(*unt);
         }
         spatial_update(unt);
      }
   }


#if TINKER_CUDART
   extern void thrust_cache_alloc();
   if (alloc_thrust_cache)
      thrust_cache_alloc();
#endif
}
TINKER_NAMESPACE_END
