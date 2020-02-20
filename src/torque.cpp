#include "torque.h"
#include "elec.h"
#include "platform.h"


TINKER_NAMESPACE_BEGIN
void torque(int vers)
{
   if (!use_elec())
      return;


#if TINKER_CUDART
      // if (pltfm_config & CU_PLTFM)
      //    torque_cu(vers);
      // else
#endif
   torque_acc(vers);
}
TINKER_NAMESPACE_END
