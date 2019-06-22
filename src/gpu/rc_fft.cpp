#include "gpu/decl_dataop.h"
#include "gpu/decl_pme.h"
#include "rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
std::vector<fft_plan_t>& fft_plans() {
  static std::vector<fft_plan_t> objs;
  return objs;
}
}
TINKER_NAMESPACE_END
