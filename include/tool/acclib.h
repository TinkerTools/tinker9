#pragma once
#include "macro.h"
#include "tool/lpflag.h"


namespace tinker {
/**
 * \ingroup rc
 * If PROCEED, directly return.
 * If DEFAULT_Q, directly return because OpenACC synchronizes the default sync
 * queue implicitly.
 * If NEW_Q, synchronize the non-default OpenACC queue with the current thread.
 */
void wait_queue(LPFlag flag);
}
