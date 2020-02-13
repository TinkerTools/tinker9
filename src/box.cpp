#include "box.h"
#include "dev_array.h"
#include "mathfunc.h"
#include "md.h"
#include "tinker_rt.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/boxes.hh>


TINKER_NAMESPACE_BEGIN
void box_data(rc_op op)
{
   if (op & rc_dealloc) {
      if (calc::traj & rc_flag) {
         box = nullptr;
         device_array::deallocate(trajbox);
      } else {
         device_array::deallocate(box);
         trajbox = nullptr;
      }
   }


   if (op & rc_alloc) {
      if (calc::traj & rc_flag) {
         device_array::allocate(trajn, &trajbox);
         box = trajbox;
      } else {
         device_array::allocate(1, &box);
      }
   }


   if (op & rc_init) {
      Box::Shape shape = Box::null;
      if (boxes::orthogonal)
         shape = Box::ortho;
      else if (boxes::monoclinic)
         shape = Box::mono;
      else if (boxes::triclinic)
         shape = Box::tri;
      else if (boxes::octahedron)
         shape = Box::oct;

      const auto& r = boxes::recip;
      const auto& l = boxes::lvec;
      recipa = make_real3(r[0][0], r[0][1], r[0][2]);
      recipb = make_real3(r[1][0], r[1][1], r[1][2]);
      recipc = make_real3(r[2][0], r[2][1], r[2][2]);
      lvec1 = make_real3(l[0][0], l[0][1], l[0][2]);
      lvec2 = make_real3(l[1][0], l[1][1], l[1][2]);
      lvec3 = make_real3(l[2][0], l[2][1], l[2][2]);

      device_array::copyin(WAIT_NEW_Q, 3, box->lvec, boxes::lvec);
      device_array::copyin(WAIT_NEW_Q, 3, box->recip, boxes::recip);
      device_array::copyin(WAIT_NEW_Q, 1, &box->volbox, &boxes::volbox);
      device_array::copyin(WAIT_NEW_Q, 1, &box->shape, &shape);
   }
}


#define DOT3(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
void copyout_box_data(const Box* pb)
{
   Box b;
   device_array::copyout(WAIT_NEW_Q, 1, &b, pb);


   if (bound::use_bounds) {
      double ax[3] = {b.lvec[0][0], b.lvec[1][0], b.lvec[2][0]};
      double bx[3] = {b.lvec[0][1], b.lvec[1][1], b.lvec[2][1]};
      double cx[3] = {b.lvec[0][2], b.lvec[1][2], b.lvec[2][2]};


      double xbox = std::sqrt(DOT3(ax, ax));
      double ybox = std::sqrt(DOT3(bx, bx));
      double zbox = std::sqrt(DOT3(cx, cx));
      double cos_a = DOT3(bx, cx) / (ybox * zbox);
      double cos_b = DOT3(cx, ax) / (zbox * xbox);
      double cos_c = DOT3(ax, bx) / (xbox * ybox);
      double a_deg = radian_dp * std::acos(cos_a);
      double b_deg = radian_dp * std::acos(cos_b);
      double c_deg = radian_dp * std::acos(cos_c);


      boxes::xbox = xbox;
      boxes::ybox = ybox;
      boxes::zbox = zbox;
      boxes::alpha = a_deg;
      boxes::beta = b_deg;
      boxes::gamma = c_deg;
      TINKER_RT(lattice)();


      const auto& r = boxes::recip;
      const auto& l = boxes::lvec;
      recipa = make_real3(r[0][0], r[0][1], r[0][2]);
      recipb = make_real3(r[1][0], r[1][1], r[1][2]);
      recipc = make_real3(r[2][0], r[2][1], r[2][2]);
      lvec1 = make_real3(l[0][0], l[0][1], l[0][2]);
      lvec2 = make_real3(l[1][0], l[1][1], l[1][2]);
      lvec3 = make_real3(l[2][0], l[2][1], l[2][2]);
   }
}


real volbox()
{
   return lvec1.x * lvec2.y * lvec3.z;
}
TINKER_NAMESPACE_END
