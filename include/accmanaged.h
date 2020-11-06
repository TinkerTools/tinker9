#pragma once
#include "mdprec.h"
#include "tool/energy_buffer.h"
#include "tool/rc_man.h"


namespace tinker {
void accmanaged_data(rc_op);


namespace detail {
extern energy_buffer_traits::type host_e_val;
extern energy_buffer_traits::type host_e_vdw;
extern energy_buffer_traits::type host_e_ele;
extern void* dptr_e_val;
extern void* dptr_e_vdw;
extern void* dptr_e_ele;
}
}
