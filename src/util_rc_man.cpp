#include "util_rc_man.h"

// GNU Fortran
#ifdef TINKER_GFORTRAN
extern "C" void _gfortran_set_args(int, char**);
#endif

// Intel
#ifdef TINKER_IFORT
extern "C" void for_rtl_init_(int*, char**);
extern "C" void for_rtl_finish_();
#endif

TINKER_NAMESPACE_BEGIN
void fortran_runtime_initialize(int argc, char** argv) {
#ifdef TINKER_GFORTRAN
  _gfortran_set_args(argc, argv);
#endif

#ifdef TINKER_IFORT
  for_rtl_init_(&argc, argv);
#endif
}

void fortran_runtime_finish() {
#ifdef TINKER_IFORT
  for_rtl_finish_();
#endif
}

void tinker_gpu_runtime_initialize() {
  const rc_t rc = static_cast<rc_t>(rc_alloc | rc_copyin);
  host_data(rc);
  device_data(rc);
}

void tinker_gpu_runtime_finish() {
  const rc_t rc = rc_dealloc;
  device_data(rc);
  host_data(rc);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
extern void random_data(rc_op);
void host_data(rc_op rc) { rc_man rand42_{random_data, rc}; }

namespace gpu {
extern void n_data(rc_t);

extern void xyz_data(rc_t);
extern void vel_data(rc_t);
extern void mass_data(rc_t);

extern void potential_data(rc_t);

extern void box_data(rc_t);
extern void couple_data(rc_t);
extern void nblist_data(rc_t);

extern void md_data(rc_t);
}

void device_data(rc_t rc) {
  using namespace gpu;

  rc_man n42_{n_data, rc};

  rc_man xyz42_{xyz_data, rc};
  rc_man vel42_{vel_data, rc};
  rc_man mass42_{mass_data, rc};

  rc_man pd42_{potential_data, rc};

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
  rc_man box42_{box_data, rc};
  rc_man cpl42_{couple_data, rc};
  rc_man nbl42_{nblist_data, rc};

  rc_man md42_{md_data, rc};
}

TINKER_NAMESPACE_END
