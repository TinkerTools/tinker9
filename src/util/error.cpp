#include "util/error.h"

TINKER_NAMESPACE_BEGIN
FatalError::FatalError(const char* _msg) : msg_(_msg) {}
FatalError::FatalError(const std::string& _msg) : msg_(_msg) {}
FatalError::FatalError(const FatalError& _e) : msg_(_e.msg_) {}
const char* FatalError::what() const noexcept { return msg_.c_str(); }
TINKER_NAMESPACE_END
