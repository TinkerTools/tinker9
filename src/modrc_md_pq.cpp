#include "mod_box.h"
#include "mod_md.h"
#include "util_array.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN

//======================================================================
// number of atoms

void n_data(rc_op op) {
  if (op & rc_dealloc) {
    trajn = -1;
    n = 0;
  }

  if (op & rc_alloc) {
    n = atoms::n;

    if (calc::traj & use_data) {
      // trajn must have been initialized by this point
      assert(trajn >= 0);
    }
  }
}

//======================================================================
// x y z coordinates

void xyz_data(rc_op op) {
  if ((calc::xyz & use_data) == 0)
    return;

  if (op & rc_dealloc) {
    if (calc::traj & use_data) {
      dealloc_array(trajx);
      dealloc_array(trajy);
      dealloc_array(trajz);
      x = nullptr;
      y = nullptr;
      z = nullptr;
    } else {
      trajx = nullptr;
      trajy = nullptr;
      trajz = nullptr;
      dealloc_array(x);
      dealloc_array(y);
      dealloc_array(z);
    }
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    if (calc::traj & use_data) {
      size *= trajn;
      alloc_array(&trajx, size);
      alloc_array(&trajy, size);
      alloc_array(&trajz, size);
      x = trajx;
      y = trajy;
      z = trajz;
    } else {
      alloc_array(&x, size);
      alloc_array(&y, size);
      alloc_array(&z, size);
    }
  }

  if (op & rc_init) {
    copyin_array(x, atoms::x, n);
    copyin_array(y, atoms::y, n);
    copyin_array(z, atoms::z, n);
  }
}

//======================================================================
// velocitiies

void vel_data(rc_op op) {
  if ((calc::vel & use_data) == 0)
    return;

  if (op & rc_dealloc) {
    dealloc_array(vx);
    dealloc_array(vy);
    dealloc_array(vz);
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    alloc_array(&vx, size);
    alloc_array(&vy, size);
    alloc_array(&vz, size);
  }

  if (op & rc_init) {
    copyin_array2(0, 3, vx, moldyn::v, n);
    copyin_array2(1, 3, vy, moldyn::v, n);
    copyin_array2(2, 3, vz, moldyn::v, n);
  }
}

//======================================================================
// atomic mass

void mass_data(rc_op op) {
  if ((calc::mass & use_data) == 0)
    return;

  if (op & rc_dealloc) {
    dealloc_array(mass);
    dealloc_array(massinv);
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    alloc_array(&mass, size);
    alloc_array(&massinv, size);
  }

  if (op & rc_init) {
    copyin_array(mass, atomid::mass, n);
    std::vector<double> mbuf(n);
    for (int i = 0; i < n; ++i)
      mbuf[i] = 1 / atomid::mass[i];
    copyin_array(massinv, mbuf.data(), n);
  }
}

void goto_frame0(int idx0) {
  assert(calc::traj & use_data);
  x = trajx + n * idx0;
  y = trajy + n * idx0;
  z = trajz + n * idx0;
  box = trajbox + idx0;
}

void goto_frame1(int idx1) { goto_frame0(idx1 - 1); }

TINKER_NAMESPACE_END
