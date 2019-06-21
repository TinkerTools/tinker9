#include "gpu/gpu.h"
#include "gpu/decl_dataop.h"

m_tinker_using_namespace;
using namespace gpu;

extern "C" {
void tinker_gpu_data_create() {
  const int op = op_alloc | op_copyin;
  gpu_data(op);
}

void tinker_gpu_data_create_() { tinker_gpu_data_create_(); }

void tinker_gpu_data_destroy() {
  const int op = op_dealloc;
  gpu_data(op);
}

void tinker_gpu_data_destroy_() { tinker_gpu_data_destroy(); }
}
