#pragma once
#include "dmflag.h"


namespace tinker {
/**
 * \brief
 * If PROCEED, directly return.
 * If DEFAULT_Q, directly return because OpenACC synchronizes the default sync
 * queue implicitly.
 * If NEW_Q, synchronize the non-default OpenACC queue with the current thread.
 */
void wait_queue(DMFlag flag);
}
