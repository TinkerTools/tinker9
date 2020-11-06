#include "accmanaged.h"
#if _OPENACC
#   include <openacc.h>
#endif


namespace tinker {
void accmanaged_data(rc_op op)
{
   using namespace detail;


   if (op & rc_dealloc) {
      dptr_e_val = nullptr;
      dptr_e_ele = nullptr;
      dptr_e_vdw = nullptr;
      #pragma acc exit data async delete(host_e_val,host_e_vdw,host_e_ele)
   }


   if (op & rc_alloc) {
      dptr_e_val = &host_e_val;
      dptr_e_vdw = &host_e_vdw;
      dptr_e_ele = &host_e_ele;
      #pragma acc enter data async create(host_e_val,host_e_vdw,host_e_ele)
#if _OPENACC
      dptr_e_val = acc_deviceptr(&host_e_val);
      dptr_e_vdw = acc_deviceptr(&host_e_vdw);
      dptr_e_ele = acc_deviceptr(&host_e_ele);
#endif
   }
}


namespace detail {
energy_buffer_traits::type host_e_val;
energy_buffer_traits::type host_e_vdw;
energy_buffer_traits::type host_e_ele;
void* dptr_e_val;
void* dptr_e_vdw;
void* dptr_e_ele;
}
}
