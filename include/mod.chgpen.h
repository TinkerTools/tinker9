#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
// enum class chpdamp_t
// {
//    gordon1,    
//    gordon2  
// };
// TINKER_EXTERN chpdamp_t pentyp;
TINKER_EXTERN real* pcore;
TINKER_EXTERN real* pval;
TINKER_EXTERN real* pval0;
TINKER_EXTERN real* palplha;
// TINKER_EXTERN count_buffer ncp;


//====================================================================//


TINKER_EXTERN real w2scale;
TINKER_EXTERN real w3scale;
TINKER_EXTERN real w4scale;
TINKER_EXTERN real w5scale;

TINKER_EXTERN int ndwexclude;
TINKER_EXTERN int (*dwexclude)[2];
TINKER_EXTERN real (*dwexclude_scale)[2];


TINKER_EXTERN int ndexclude;
TINKER_EXTERN int (*dexclude)[2];
TINKER_EXTERN real* dexclude_scale;


TINKER_EXTERN int nwexclude;
TINKER_EXTERN int (*wexclude)[2];
TINKER_EXTERN real* wexclude_scale;

}
