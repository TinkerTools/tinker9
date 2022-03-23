#include "tool/rcman.h"

namespace tinker {
bool ResourceManagement::will_dealloc() const
{
   return m_op & rc_dealloc;
}

bool ResourceManagement::only_dealloc() const
{
   return m_op == rc_dealloc;
}

ResourceManagement::ResourceManagement(void (*f)(RcOp), RcOp op)
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
   RcOp op = rc_alloc | rc_init;
   device_data(op);
}

void finish()
{
   RcOp op = rc_dealloc;
   device_data(op);
}
}

#include "math/inc.h"
#include "platform.h"
#include "tool/gpucard.h"

#include "ff/box.h"
#include "ff/energy.h"
#include "ff/molecule.h"
#include "ff/nblist.h"
#include "ff/rattle.h"
#include "md/inc.h"
#include "md/osrw.h"
#include "tool/cudalib.h"

namespace tinker {
void device_data(RcOp op)
{
   // host
   RcMan rand42{randomData, op};
   RcMan pf42{platformData, op};
   RcMan gpu42{gpuData, op};

   // device
   RcMan cl42{cudalibData, op};

   RcMan n42{mdNData, op};
   RcMan cpl42{coupleData, op};

   RcMan box42{boxData, op};

   RcMan xyz42{mdXyzData, op};
   RcMan vel42{mdVelData, op};
   RcMan mass42{mdMassData, op};
   RcMan molecule42{moleculeData, op};
   RcMan group42{groupData, op};

   RcMan energy42{energyData, op};
   RcMan osrw42{osrw_data, op};

   // Neighbor lists must be initialized after potential initialization.
   // xred, yred, and zred need to be initialized in vdw (Halgren 14-7)
   // and will be used in nblist setups.
   RcMan nbl42{nblist_data, op};

   RcMan rattle42{rattleData, op};
   RcMan md42{mdData, op};
}
}
