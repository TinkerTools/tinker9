#include "array.h"
#include "box.h"
#include "md.h"
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

    if (calc::traj & rc_flag) {
      // trajn must have been initialized by this point
      assert(trajn >= 0);
    }
  }
}

//======================================================================
// x y z coordinates

void xyz_data(rc_op op) {
  if ((calc::xyz & rc_flag) == 0)
    return;

  if (op & rc_dealloc) {
    if (calc::traj & rc_flag) {
      dealloc_bytes(trajx);
      dealloc_bytes(trajy);
      dealloc_bytes(trajz);
      x = nullptr;
      y = nullptr;
      z = nullptr;
    } else {
      trajx = nullptr;
      trajy = nullptr;
      trajz = nullptr;
      dealloc_bytes(x);
      dealloc_bytes(y);
      dealloc_bytes(z);
    }
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    if (calc::traj & rc_flag) {
      size *= trajn;
      alloc_bytes(&trajx, size);
      alloc_bytes(&trajy, size);
      alloc_bytes(&trajz, size);
      x = trajx;
      y = trajy;
      z = trajz;
    } else {
      alloc_bytes(&x, size);
      alloc_bytes(&y, size);
      alloc_bytes(&z, size);
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
  if ((calc::vel & rc_flag) == 0)
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(vx);
    dealloc_bytes(vy);
    dealloc_bytes(vz);
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    alloc_bytes(&vx, size);
    alloc_bytes(&vy, size);
    alloc_bytes(&vz, size);
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
  if ((calc::mass & rc_flag) == 0)
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(mass);
    dealloc_bytes(massinv);
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    alloc_bytes(&mass, size);
    alloc_bytes(&massinv, size);
  }

  if (op & rc_init) {
    copyin_array(mass, atomid::mass, n);
    std::vector<double> mbuf(n);
    for (int i = 0; i < n; ++i)
      mbuf[i] = 1 / atomid::mass[i];
    copyin_array(massinv, mbuf.data(), n);
  }
}

void goto_frame(int idx0) {
  assert(calc::traj & rc_flag);
  x = trajx + n * idx0;
  y = trajy + n * idx0;
  z = trajz + n * idx0;
  box = trajbox + idx0;
}
TINKER_NAMESPACE_END
