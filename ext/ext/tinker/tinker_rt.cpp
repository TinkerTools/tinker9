#include "tinker_rt.h"
#include "tinker_mod.h"

TINKER_NAMESPACE_BEGIN
void nextarg(size_t len, char* str, logical& exist) {
  const char blank = ' ';
  std::memset(str, blank, len);
  exist = _false_;

  if (argue::narg != 0) {
    size_t length = std::min(len, sizeof(argue::arg[argue::maxarg]));
    for (size_t i = 1; i <= argue::narg; ++i) {
      if (argue::listarg[i]) {
        argue::listarg[i] = _false_;
        std::strncpy(str, argue::arg[i], length);
        exist = _true_;
        break;
      }
    }
  }
}
TINKER_NAMESPACE_END
