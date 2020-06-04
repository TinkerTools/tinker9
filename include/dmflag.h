#pragma once
#include "enum_op.h"


namespace tinker {
/**
 * \ingroup rc
 * Policy of launching kernels.
 */
enum class DMFlag
{
   DEFAULT_Q = 0x01, ///< Default CUDA stream / OpenACC queue.
   NEW_Q = 0x00,     ///< Non-default CUDA stream / OpenACC queue.
   WAIT = 0x02,      ///< Synchronous / Blocking.
   PROCEED = 0x00,   ///< Asynchronous / Non-Blocking.
};
TINKER_ENABLE_ENUM_BITMASK(DMFlag);


/**
 * \ingroup rc
 * Executes kernel on the default CUDA stream / OpenACC queue. Non-blocking.
 */
constexpr DMFlag PROCEED_DEFAULT_Q = (DMFlag::PROCEED | DMFlag::DEFAULT_Q);
/**
 * \ingroup rc
 * Executes kernel on the non-default CUDA stream / OpenACC queue. Non-blocking.
 */
constexpr DMFlag PROCEED_NEW_Q = (DMFlag::PROCEED | DMFlag::NEW_Q);
/**
 * \ingroup rc
 * Executes kernel on the default CUDA stream / OpenACC queue. Blocking.
 */
constexpr DMFlag WAIT_DEFAULT_Q = (DMFlag::WAIT | DMFlag::DEFAULT_Q);
/**
 * \ingroup rc
 * Executes kernel on the non-default CUDA stream / OpenACC queue. Blocking.
 */
constexpr DMFlag WAIT_NEW_Q = (DMFlag::WAIT | DMFlag::NEW_Q);
}
