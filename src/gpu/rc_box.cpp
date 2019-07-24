#include "gpu/decl_box.h"
#include "gpu/decl_mdstate.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void box_data(rc_t rc) {
  if (rc & rc_dealloc) {
    if (use_traj & use_data) {
      box = nullptr;
      check_cudart(cudaFree(trajbox));
    } else {
      check_cudart(cudaFree(box));
      trajbox = nullptr;
    }
  }

  if (rc & rc_alloc) {
    size_t size = sizeof(box_t);
    if (use_traj & use_data) {
      size *= trajn;
      check_cudart(cudaMalloc(&trajbox, size));
      box = trajbox;
    } else {
      check_cudart(cudaMalloc(&box, size));
    }
  }

  if (rc & rc_copyin) {
    box_t::shape_t shape = box_t::null;
    if (boxes::orthogonal)
      shape = box_t::ortho;
    else if (boxes::monoclinic)
      shape = box_t::mono;
    else if (boxes::triclinic)
      shape = box_t::tri;
    else if (boxes::octahedron)
      shape = box_t::oct;

    copyin_array(&box->lvec[0][0], &boxes::lvec[0][0], 9);
    copyin_array(&box->recip[0][0], &boxes::recip[0][0], 9);
    copyin_array(&box->volbox, &boxes::volbox, 1);
    check_cudart(cudaMemcpy(&box->shape, &shape, sizeof(box_t::shape_t),
                            cudaMemcpyHostToDevice));
  }
}

void box_data_copyout(const box_t& b) {
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
}
TINKER_NAMESPACE_END
