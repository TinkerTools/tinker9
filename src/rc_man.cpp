#include "rc_man.h"
#include "wait_queue.h"


namespace tinker {
bool ResourceManagement::will_dealloc() const
{
   return op_ & rc_dealloc;
}


bool ResourceManagement::only_dealloc() const
{
   return op_ == rc_dealloc;
}


ResourceManagement::ResourceManagement(void (*f)(rc_op), rc_op op)
   : f_(f)
   , op_(op)
{
   if (!will_dealloc()) {
      f_(op_);
   }
}


ResourceManagement::~ResourceManagement()
{
   if (only_dealloc()) {
      f_(op_);
   }
}


void initialize()
{
   rc_op op = rc_alloc | rc_init;
   device_data(op);
}


void finish()
{
   rc_op op = rc_dealloc;
   device_data(op);
}
}


#include "gpu_card.h"
#include "platform.h"
#include "random.h"


#include "box.h"
#include "couple.h"
#include "cudalib.h"
#include "energy.h"
#include "group.h"
#include "mdintg.h"
#include "mdpq.h"
#include "molecule.h"
#include "nblist.h"
#include "osrw.h"


namespace tinker {
void device_data(rc_op op)
{
   // host
   rc_man rand42{random_data, op};
   rc_man pf42{platform_data, op};
   rc_man gpu_card42{gpu_card_data, op};


   // device
   rc_man cl42{cudalib_data, op};

   rc_man n42{n_data, op};
   rc_man cpl42{couple_data, op};

   rc_man box42{box_data, op};

   rc_man xyz42{xyz_data, op};
   rc_man vel42{vel_data, op};
   rc_man mass42{mass_data, op};
   rc_man molecule42{molecule_data, op};
   rc_man group42{group_data, op};

   rc_man energy42{energy_data, op};
   rc_man osrw42{osrw_data, op};

   // Neighbor lists must be initialized after potential initialization.
   // xred, yred, and zred need to be initialized in vdw (Halgren 14-7)
   // and will be used in nblist setups.
   rc_man nbl42{nblist_data, op};

   rc_man md42{md_data, op};
}
}


#if defined(TINKER_GFORTRAN)
// GNU Fortran
extern "C" void _gfortran_set_args(int, char**);


namespace tinker {
void fortran_runtime_initialize(int argc, char** argv)
{
   _gfortran_set_args(argc, argv);
}


void fortran_runtime_finish() {}
}


#elif defined(TINKER_IFORT)
// Intel
extern "C" void for_rtl_init_(int*, char**);
extern "C" int for_rtl_finish_();
extern "C" int for__l_argc;
extern "C" char** for__a_argv;


namespace tinker {
void fortran_runtime_initialize(int argc, char** argv)
{
   // what a surprise that I had to use these undocumented variables here
   for__l_argc = argc;
   for__a_argv = argv;
   for_rtl_init_(&argc, argv);
}


void fortran_runtime_finish()
{
   for_rtl_finish_();
}
}
#endif
