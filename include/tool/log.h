#pragma once
#include "tool/error.h"


#ifndef TINKER_ENABLE_LOG
#   define TINKER_ENABLE_LOG 0
#endif
#if TINKER_ENABLE_LOG
#   define TINKER_LOG(...)                                                     \
      print(stderr, " LOG -- %s (at %s:%d)\n", format(__VA_ARGS__), __FILE__,  \
            __LINE__)
#else
#   define TINKER_LOG(...)
#endif
