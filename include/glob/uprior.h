#pragma once
#include "macro.h"

// uprior
namespace tinker {
enum class UPred
{
   NONE,
   GEAR,
   ASPC,
   LSQR
};
TINKER_EXTERN UPred polpred;
TINKER_EXTERN int maxualt;
TINKER_EXTERN int nualt;
TINKER_EXTERN real (*udalt_00)[3];
TINKER_EXTERN real (*udalt_01)[3];
TINKER_EXTERN real (*udalt_02)[3];
TINKER_EXTERN real (*udalt_03)[3];
TINKER_EXTERN real (*udalt_04)[3];
TINKER_EXTERN real (*udalt_05)[3];
TINKER_EXTERN real (*udalt_06)[3];
TINKER_EXTERN real (*udalt_07)[3];
TINKER_EXTERN real (*udalt_08)[3];
TINKER_EXTERN real (*udalt_09)[3];
TINKER_EXTERN real (*udalt_10)[3];
TINKER_EXTERN real (*udalt_11)[3];
TINKER_EXTERN real (*udalt_12)[3];
TINKER_EXTERN real (*udalt_13)[3];
TINKER_EXTERN real (*udalt_14)[3];
TINKER_EXTERN real (*udalt_15)[3];
TINKER_EXTERN real (*upalt_00)[3];
TINKER_EXTERN real (*upalt_01)[3];
TINKER_EXTERN real (*upalt_02)[3];
TINKER_EXTERN real (*upalt_03)[3];
TINKER_EXTERN real (*upalt_04)[3];
TINKER_EXTERN real (*upalt_05)[3];
TINKER_EXTERN real (*upalt_06)[3];
TINKER_EXTERN real (*upalt_07)[3];
TINKER_EXTERN real (*upalt_08)[3];
TINKER_EXTERN real (*upalt_09)[3];
TINKER_EXTERN real (*upalt_10)[3];
TINKER_EXTERN real (*upalt_11)[3];
TINKER_EXTERN real (*upalt_12)[3];
TINKER_EXTERN real (*upalt_13)[3];
TINKER_EXTERN real (*upalt_14)[3];
TINKER_EXTERN real (*upalt_15)[3];

TINKER_EXTERN real* udalt_lsqr_a;
TINKER_EXTERN real* upalt_lsqr_a;
TINKER_EXTERN real* udalt_lsqr_b;
TINKER_EXTERN real* upalt_lsqr_b;
}
