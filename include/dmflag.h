#pragma once
#include "enum_op.h"


TINKER_NAMESPACE_BEGIN
enum class DMFlag
{
   DEFAULT_Q = 0x01, // vs. NEW_Q
   NEW_Q = 0x00,

   WAIT = 0x02, // vs. PROCEED
   PROCEED = 0x00,
};
TINKER_ENABLE_ENUM_BITMASK(DMFlag);


constexpr DMFlag PROCEED_DEFAULT_Q = (DMFlag::PROCEED | DMFlag::DEFAULT_Q);
constexpr DMFlag PROCEED_NEW_Q = (DMFlag::PROCEED | DMFlag::NEW_Q);
constexpr DMFlag WAIT_DEFAULT_Q = (DMFlag::WAIT | DMFlag::DEFAULT_Q);
constexpr DMFlag WAIT_NEW_Q = (DMFlag::WAIT | DMFlag::NEW_Q);
TINKER_NAMESPACE_END
