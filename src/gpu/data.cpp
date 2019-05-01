#include "gpu/data.h"
#include "gpu/mdstate.h"

extern "C" {
void tinker_gpu_data_create() {
  m_tinker_using_namespace;
  using namespace gpu;
  accel_data(op_create);
  mass_data(op_create);
  energy_data(op_create);
}

void tinker_gpu_data_destroy() {
  m_tinker_using_namespace;
  using namespace gpu;
  xyz_data(op_destroy);
  vel_data(op_destroy);
  accel_data(op_destroy);
  mass_data(op_destroy);
  energy_data(op_destroy);
}
}
