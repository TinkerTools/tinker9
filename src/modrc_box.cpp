#include "mod_box.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_math.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

TINKER_NAMESPACE_BEGIN
void box_data(rc_op op) {
  if (op & rc_dealloc) {
    if (calc::traj & use_data) {
      box = nullptr;
      dealloc_array(trajbox);
    } else {
      dealloc_array(box);
      trajbox = nullptr;
    }
  }

  if (op & rc_alloc) {
    size_t size = sizeof(Box);
    if (calc::traj & use_data) {
      size *= trajn;
      alloc_array(&trajbox, size);
      box = trajbox;
    } else {
      alloc_array(&box, size);
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

    copyin_array(&box->lvec[0][0], &boxes::lvec[0][0], 9);
    copyin_array(&box->recip[0][0], &boxes::recip[0][0], 9);
    copyin_array(&box->volbox, &boxes::volbox, 1);
    copy_memory(&box->shape, &shape, sizeof(Box::Shape),
                CopyDirection::HostToDevice);
  }
}

void box_data_copyout(const Box& b) {
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
