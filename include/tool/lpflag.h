#pragma once
#include "tool/enum_op.h"


namespace tinker {
/**
 * \ingroup rc
 * Launching policies.
 */
enum class LPFlag
{
   DEFAULT_Q = 0x01, ///< Default CUDA stream / OpenACC queue.
   NEW_Q = 0x02,     ///< Non-default CUDA stream / OpenACC queue.
   WAIT = 0x04,      ///< Synchronous / Blocking.
   PROCEED = 0x08,   ///< Asynchronous / Non-Blocking.
};
TINKER_ENABLE_ENUM_BITMASK(LPFlag);


/**
 * \ingroup rc
 * Executes kernel on the default CUDA stream / OpenACC queue. Non-blocking.
 */
constexpr LPFlag PROCEED_DEFAULT_Q = (LPFlag::PROCEED | LPFlag::DEFAULT_Q);
/**
 * \ingroup rc
 * Executes kernel on the non-default CUDA stream / OpenACC queue. Non-blocking.
 */
constexpr LPFlag PROCEED_NEW_Q = (LPFlag::PROCEED | LPFlag::NEW_Q);
/**
 * \ingroup rc
 * Executes kernel on the default CUDA stream / OpenACC queue. Blocking.
 */
constexpr LPFlag WAIT_DEFAULT_Q = (LPFlag::WAIT | LPFlag::DEFAULT_Q);
/**
 * \ingroup rc
 * Executes kernel on the non-default CUDA stream / OpenACC queue. Blocking.
 */
constexpr LPFlag WAIT_NEW_Q = (LPFlag::WAIT | LPFlag::NEW_Q);
}
