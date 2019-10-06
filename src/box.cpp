#include "box.h"
#include "dev_array.h"
#include "mathfunc.h"
#include "md.h"
#include "tinker_rt.h"
#include <ext/tinker/detail/bound.hh>
#include <ext/tinker/detail/boxes.hh>

TINKER_NAMESPACE_BEGIN
void box_data(rc_op op) {
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

    device_array::copyin(3, box->lvec, boxes::lvec);
    device_array::copyin(3, box->recip, boxes::recip);
    device_array::copyin(1, &box->volbox, &boxes::volbox);
    device_array::copyin(1, &box->shape, &shape);
  }
}

void copyout_box_data(const Box* pb) {
  Box b;
  device_array::copyout(1, &b, pb);

  if (bound::use_bounds) {
    double ax[3] = {b.lvec[0][0], b.lvec[1][0], b.lvec[2][0]};
    double bx[3] = {b.lvec[0][1], b.lvec[1][1], b.lvec[2][1]};
    double cx[3] = {b.lvec[0][2], b.lvec[1][2], b.lvec[2][2]};

#define DOT_IMPL_(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
    double xbox = std::sqrt(DOT_IMPL_(ax, ax));
    double ybox = std::sqrt(DOT_IMPL_(bx, bx));
    double zbox = std::sqrt(DOT_IMPL_(cx, cx));
    double cos_a = DOT_IMPL_(bx, cx) / (ybox * zbox);
    double cos_b = DOT_IMPL_(cx, ax) / (zbox * xbox);
    double cos_c = DOT_IMPL_(ax, bx) / (xbox * ybox);
    double a_deg = radian * std::acos(cos_a);
    double b_deg = radian * std::acos(cos_b);
    double c_deg = radian * std::acos(cos_c);
#undef DOT_IMPL_

    boxes::xbox = xbox;
    boxes::ybox = ybox;
    boxes::zbox = zbox;
    boxes::alpha = a_deg;
    boxes::beta = b_deg;
    boxes::gamma = c_deg;
    TINKER_RT(lattice)();
  }
}
TINKER_NAMESPACE_END
