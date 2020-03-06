#include "rc_man.h"
#include "wait_queue.h"


TINKER_NAMESPACE_BEGIN
bool ResourceManagement::will_dealloc_() const
{
   return op_ & rc_dealloc;
}


bool ResourceManagement::only_dealloc_() const
{
   return op_ == rc_dealloc;
}


ResourceManagement::ResourceManagement(void (*f)(rc_op), rc_op op)
   : f_(f)
   , op_(op)
{
   if (!will_dealloc_()) {
      f_(op_);
   }
}


ResourceManagement::~ResourceManagement()
{
   if (only_dealloc_()) {
      f_(op_);
   }
}


void initialize()
{
   rc_op op = rc_alloc | rc_init;
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


#include "gpu_card.h"
#include "platform.h"
#include "random.h"


#include "box.h"
#include "cudalib.h"
#include "energy.h"
#include "group.h"
#include "mdintg.h"
#include "mdpq.h"
#include "molecule.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
void host_data(rc_op op)
{
   rc_man rand42_{random_data, op};
   rc_man pf42_{platform_data, op};
   rc_man gpu_card42_{gpu_card_data, op};
}


void device_data(rc_op op)
{
   rc_man cl42_{cudalib_data, op};

   rc_man n42_{n_data, op};

   rc_man box42_{box_data, op};

   rc_man xyz42_{xyz_data, op};
   rc_man vel42_{vel_data, op};
   rc_man mass42_{mass_data, op};
   rc_man molecule42_{molecule_data, op};
   rc_man group42_{group_data, op};

   rc_man energy42_{energy_data, op};

   // Neighbor lists must be initialized after potential initialization.
   // xred, yred, and zred need to be initialized in vdw (Halgren 14-7)
   // and will be used in nblist setups.
   rc_man nbl42_{nblist_data, op};

   rc_man md42_{md_data, op};
}
TINKER_NAMESPACE_END


#if defined(TINKER_GFORTRAN)
// GNU Fortran
extern "C" void _gfortran_set_args(int, char**);


TINKER_NAMESPACE_BEGIN
void fortran_runtime_initialize(int argc, char** argv)
{

   _gfortran_set_args(argc, argv);
}


void fortran_runtime_finish() {}
TINKER_NAMESPACE_END


#elif defined(TINKER_IFORT)
// Intel
extern "C" void for_rtl_init_(int*, char**);
extern "C" void for_rtl_finish_();


TINKER_NAMESPACE_BEGIN
void fortran_runtime_initialize(int argc, char** argv)
{
   for_rtl_init_(&argc, argv);
}


void fortran_runtime_finish()
{
   for_rtl_finish_();
}
TINKER_NAMESPACE_END
#endif
