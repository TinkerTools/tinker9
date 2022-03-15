#include "tool/rc_man.h"

namespace tinker {
bool ResourceManagement::will_dealloc() const
{
   return m_op & rc_dealloc;
}

bool ResourceManagement::only_dealloc() const
{
   return m_op == rc_dealloc;
}

ResourceManagement::ResourceManagement(void (*f)(rc_op), rc_op op)
   : m_f(f)
   , m_op(op)
{
   if (!will_dealloc()) {
      m_f(m_op);
   }
}

ResourceManagement::~ResourceManagement()
{
   if (only_dealloc()) {
      m_f(m_op);
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

#include "platform.h"
#include "random.h"
#include "tool/gpu_card.h"

#include "box.h"
#include "couple.h"
#include "energy.h"
#include "group.h"
#include "mdintg.h"
#include "mdpq.h"
#include "molecule.h"
#include "nblist.h"
#include "osrw.h"
#include "rattle.h"
#include "tool/cudalib.h"

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

   rc_man rattle42{rattle_data, op};
   rc_man md42{md_data, op};
}
}
