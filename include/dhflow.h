#pragma once
#include "tool/energy_buffer.h"


namespace tinker {
struct DHFlow
{
   energy_buffer_traits::type e_val;
   energy_buffer_traits::type e_vdw;
   energy_buffer_traits::type e_ele;
   virial_buffer_traits::type v_val[virial_buffer_traits::N];
   virial_buffer_traits::type v_vdw[virial_buffer_traits::N];
   virial_buffer_traits::type v_ele[virial_buffer_traits::N];
};
}
