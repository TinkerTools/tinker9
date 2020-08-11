#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
enum class chpdamp_t
{
   gordon1,    
   gordon2  
};
TINKER_EXTERN chpdamp_t pentyp;
TINKER_EXTERN real* pcore;
TINKER_EXTERN real* pval;
TINKER_EXTERN real* pval0;
TINKER_EXTERN real* palplha;
TINKER_EXTERN count_buffer ncp;
}
