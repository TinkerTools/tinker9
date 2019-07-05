#include "util/error.h"

TINKER_NAMESPACE_BEGIN
FatalError::FatalError(const char* msg) : msg_(msg) {}
FatalError::FatalError(const std::string& msg) : msg_(msg) {}
FatalError::FatalError(const FatalError& e) : msg_(e.msg_) {}
const char* FatalError::what() const noexcept { return msg_.c_str(); }
TINKER_NAMESPACE_END
