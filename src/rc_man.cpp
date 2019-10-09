#include "rc_man.h"

#ifdef TINKER_GFORTRAN
// GNU Fortran
extern "C" void _gfortran_set_args(int, char**);
#elif TINKER_IFORT
// Intel
extern "C" void for_rtl_init_(int*, char**);
extern "C" void for_rtl_finish_();
#else
#   error "unknown fortran compiler error"
#   error see also "macro.h"
#endif

TINKER_NAMESPACE_BEGIN
void fortran_runtime_initialize(int argc, char** argv)
{
#ifdef TINKER_GFORTRAN
   _gfortran_set_args(argc, argv);
#elif TINKER_IFORT
   for_rtl_init_(&argc, argv);
#else
#   error "unknown fortran compiler error"
#   error see also "macro.h"
#endif
}

void fortran_runtime_finish()
{
#ifdef TINKER_GFORTRAN
#elif TINKER_IFORT
   for_rtl_finish_();
#else
#   error "unknown fortran compiler error"
#   error see also "macro.h"
#endif
}

void initialize()
{
   rc_op op = static_cast<rc_op>(rc_alloc | rc_init);
   host_data(op);
   device_data(op);
}

void finish()
{
   rc_op op = rc_dealloc;
   device_data(op);
   host_data(op);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
extern void random_data(rc_op);
extern void gpu_card_data(rc_op);
void host_data(rc_op op)
{
   rc_man rand42_{random_data, op};
   rc_man gpu_card42_{gpu_card_data, op};
}

extern void n_data(rc_op);

extern void xyz_data(rc_op);
extern void vel_data(rc_op);
extern void mass_data(rc_op);

extern void potential_data(rc_op);

extern void box_data(rc_op);
extern void nblist_data(rc_op);

extern void md_data(rc_op);

void device_data(rc_op op)
{
   rc_man n42_{n_data, op};

   rc_man xyz42_{xyz_data, op};
   rc_man vel42_{vel_data, op};
   rc_man mass42_{mass_data, op};

   rc_man pd42_{potential_data, op};

   // Neighbor lists must be initialized after potential initialization.
   // xred, yred, and zred need to be initialized in vdw routines and will be
   // used in nblist setups.
   rc_man box42_{box_data, op};
   rc_man nbl42_{nblist_data, op};

   rc_man md42_{md_data, op};
}

TINKER_NAMESPACE_END
