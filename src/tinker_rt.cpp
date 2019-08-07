#include "tinker_rt.h"
#include <algorithm>
#include <ext/tinker/detail/argue.hh>

TINKER_NAMESPACE_BEGIN
void nextarg(size_t len, char* str, int& exist) {
  const char blank = ' ';
  std::memset(str, blank, len);
  exist = false;

  if (argue::narg != 0) {
    size_t length = std::min(len, sizeof(argue::arg[argue::maxarg]));
    for (size_t i = 1; i <= argue::narg; ++i) {
      if (argue::listarg[i]) {
        argue::listarg[i] = false;
        std::strncpy(str, argue::arg[i], length);
        exist = true;
        break;
      }
    }
  }
}
TINKER_NAMESPACE_END
