#pragma once
#include "glob.gpucard.h"
#include "tool/rc_man.h"


namespace tinker {
void gpu_card_data(rc_op op);
int get_grid_size(int nthreads_per_block);
}


#if TINKER_CUDART
namespace tinker {
std::string get_cuda_runtime_version_string();
std::string get_cuda_driver_version_string();
std::string get_thrust_version_string();
std::vector<DeviceAttribute>& get_device_attributes();
}
#endif
