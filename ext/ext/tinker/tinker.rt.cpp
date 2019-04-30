#include "tinker.rt.h"
#include "tinker.mod.h"

TINKER_NAMESPACE_BEGIN
void nextarg(size_t _len, char* _str, logical& _exist) {
  const char blank = ' ';
  std::memset(_str, blank, _len);
  _exist = _false_;

  if (argue::narg != 0) {
    size_t length = std::min(_len, sizeof(argue::arg[argue::maxarg]));
    for (size_t i = 1; i <= argue::narg; ++i) {
      if (argue::listarg[i]) {
        argue::listarg[i] = _false_;
        std::strncpy(_str, argue::arg[i], length);
        _exist = _true_;
        break;
      }
    }
  }
}
TINKER_NAMESPACE_END
