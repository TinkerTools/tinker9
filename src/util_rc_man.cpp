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
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
extern void random_data(rc_op);
void host_data(rc_op rc) { rc_man rand42_{random_data, rc}; }
TINKER_NAMESPACE_END
