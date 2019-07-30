#include "util_pme.h"

TINKER_NAMESPACE_BEGIN
std::vector<fft_plan_t>& fft_plans() {
  static std::vector<fft_plan_t> objs;
  return objs;
}
TINKER_NAMESPACE_END
