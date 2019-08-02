#include "mod_pme.h"
#include "util_rt.h"

TINKER_NAMESPACE_BEGIN
std::vector<FFTPlan>& fft_plans() {
  static std::vector<FFTPlan> objs;
  return objs;
}
TINKER_NAMESPACE_END
