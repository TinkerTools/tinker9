#include "gpu/gpu.h"
#include "gpu/decl_mdstate.h"

m_tinker_using_namespace;
using namespace gpu;

extern "C" {
void tinker_gpu_data_create() {
  const rc_t rc = static_cast<rc_t>(rc_alloc | rc_copyin);
  mdstate_data(rc);
}

void tinker_gpu_data_create_() { tinker_gpu_data_create_(); }

void tinker_gpu_data_destroy() {
  const rc_t rc = rc_dealloc;
  mdstate_data(rc);
}

void tinker_gpu_data_destroy_() { tinker_gpu_data_destroy(); }
}
