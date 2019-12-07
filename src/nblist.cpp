#include "nblist.h"
#include "dev_array.h"
#include "e_polar.h"
#include "e_vdw.h"
#include "md.h"
#include "potent.h"
#include "spatial.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/neigh.hh>
#include <tinker/detail/potent.hh>


TINKER_NAMESPACE_BEGIN
void check_nblist_acc(int, real, const Box*, int*, const real*, const real*,
                      const real*, real*, real*, real*);
int check_spatial_acc(int, real, const Box*, int*, const real*, const real*,
                      const real*, real*, real*, real*);
int always_use_nblist;


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


// see also cutoffs.f
// In the gas phase calculation where neighbor list is not used, we should
// always first check the value of maxn.
// If maxn is equal to 1, it means the value of cutoff can even be INF.
static int nblist_maxlst_(int maxn, double cutoff, double buffer)
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

NBList::~NBList()
{
   device_array::deallocate(nlst, lst, update, xold, yold, zold);
}

static void nblist_op_alloc_(NBListUnit& nblu, int maxn, double cutoff,
                             double buffer, const real* _x, const real* _y,
                             const real* _z)
{
   nblu = NBListUnit::open();
   auto& st = *nblu;

   device_array::allocate(n, &st.nlst);

   int maxlst = nblist_maxlst_(maxn, cutoff, buffer);
   device_array::allocate(maxlst * n, &st.lst);

   if (maxlst == 1) {
      st.update = nullptr;
      st.xold = nullptr;
      st.yold = nullptr;
      st.zold = nullptr;
   } else {
      device_array::allocate(n, &st.update, &st.xold, &st.yold, &st.zold);
   }

   st.x = _x;
   st.y = _y;
   st.z = _z;

   st.maxnlst = maxlst;
   st.cutoff = cutoff;
   st.buffer = buffer;

   nblu.update_deviceptr(st);
}


extern void nblist_build_acc(NBListUnit);
extern void nblist_update_acc(NBListUnit);
extern void spatial_data_init_cu(SpatialUnit);


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


   int maxnlst = 0;
   int u = 0;
   bool alloc_thrust_cache = false;
   double cut = 0;
   double buf = 0;


   // vlist
   u = vlist_version();
   int& vunt = nblist_get_unit(u, (int&)vlist_unit, (int&)vspatial_unit);
   cut = limits::vdwcut;
   buf = neigh::lbuffer;
   if (u && (op & rc_alloc)) {
      nblist_data_alloc(u, vunt, 2500, cut, buf, xred, yred, zred);
      if (u & NBList::spatial)
         alloc_thrust_cache = true;
   }
   if (u && (op & rc_init)) {
      evdw_reduce_xyz();
      nblist_data_init(u, vunt);
   }
   if (u && (op & rc_evolve)) {
      evdw_reduce_xyz();
      nblist_data_evolve(u, vunt);
   }

   // dlist

   // clist

   // mlist
   u = mlist_version();
   int& munt = nblist_get_unit(u, (int&)mlist_unit, (int&)mspatial_unit);
   cut = limits::mpolecut;
   buf = neigh::lbuffer;
   if (u && (op & rc_alloc)) {
      if (u == NBList::spatial)
         nblist_data_alloc(NBList::nblist, (int&)mlist_unit, 2500, cut, buf, x,
                           y, z);
      nblist_data_alloc(u, munt, 2500, cut, buf, x, y, z);
      if (u & NBList::spatial)
         alloc_thrust_cache = true;
   }
   if (u && (op & rc_init)) {
      if (u == NBList::spatial)
         nblist_data_init(NBList::nblist, (int&)mlist_unit);
      nblist_data_init(u, munt);
   }
   if (u && (op & rc_evolve)) {
      if (rc_flag & calc::traj) {
         mlist_unit->x = x;
         mlist_unit->y = y;
         mlist_unit->z = z;
         mlist_unit.update_deviceptr(*mlist_unit);

         if (u == NBList::spatial) {
            mspatial_unit->x = x;
            mspatial_unit->y = y;
            mspatial_unit->z = z;
            mspatial_unit.update_deviceptr(*mspatial_unit);
         }
      }

      if (u == NBList::spatial)
         nblist_data_evolve(NBList::nblist, (int&)mlist_unit);
      nblist_data_evolve(u, munt);
   }

   // ulist
   u = ulist_version();
   if (u) {
      if (op & rc_dealloc) {
         device_array::deallocate(mindex, minv, minv_exclude_);
      }

      if (op & rc_alloc) {
         maxnlst = 500;
         int minv_size = maxnlst;
         if (u == NBList::double_loop)
            maxnlst = 1;
         nblist_op_alloc_(ulist_unit, maxnlst, limits::usolvcut, neigh::pbuffer,
                          x, y, z);
         device_array::allocate(n, &mindex);
         minv_size =
            nblist_maxlst_(minv_size, limits::usolvcut, neigh::pbuffer);
         device_array::allocate(3 * minv_size * n, &minv);
         device_array::allocate(6 * nuexclude_, &minv_exclude_);
      }

      if (op & rc_init)
         nblist_build_acc(ulist_unit);

      if (op & rc_evolve) {
         if (rc_flag & calc::traj) {
            ulist_unit->x = x;
            ulist_unit->y = y;
            ulist_unit->z = z;
            ulist_unit.update_deviceptr(*mlist_unit);
         }
         nblist_update_acc(ulist_unit);
      }
   }


#if TINKER_CUDART
   extern void thrust_cache_alloc();
   if (alloc_thrust_cache)
      thrust_cache_alloc();
#endif
}


int& nblist_get_unit(int version, int& nb, int& sp)
{
   if (version & (NBList::nblist | NBList::double_loop)) {
      return nb;
   }
#if TINKER_CUDART
   else if (version & NBList::spatial) {
      return sp;
   }
#endif
}


void nblist_data_alloc(int version, int& unt, int maxn, real cut, real buf,
                       const real* x, const real* y, const real* z)
{
   if (version & (NBList::nblist | NBList::double_loop)) {
      int maxnlst = maxn;
      if (version == NBList::double_loop)
         maxnlst = 1;
      nblist_op_alloc_((NBListUnit&)unt, maxnlst, cut, buf, x, y, z);
   }
#if TINKER_CUDART
   else if (version & NBList::spatial) {
      spatial_data_alloc((SpatialUnit&)unt, n, cut, buf, x, y, z);
   }
#endif
}


void nblist_data_init(int version, int unt)
{
   if (version & (NBList::nblist | NBList::double_loop)) {
      nblist_build_acc((NBListUnit)unt);
   }
#if TINKER_CUDART
   else if (version & NBList::spatial) {
      spatial_data_init_cu((SpatialUnit)unt);
   }
#endif
}


void nblist_data_evolve(int version, int unt)
{
   if (version & (NBList::nblist | NBList::double_loop)) {
      nblist_update_acc((NBListUnit)unt);
   }
#if TINKER_CUDART
   else if (version & NBList::spatial) {
      SpatialUnit sunt = unt;
      auto& st = *sunt;
      st.rebuild = check_spatial_acc(st.n, st.buffer, box, st.update, st.x,
                                     st.y, st.z, st.xold, st.yold, st.zold);
      if (st.rebuild)
         spatial_data_init_cu((SpatialUnit)unt);
   }
#endif
}
TINKER_NAMESPACE_END
