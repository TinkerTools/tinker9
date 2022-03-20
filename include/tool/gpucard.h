#pragma once
#include "tool/rcman.h"

namespace tinker {
void gpuData(RcOp);
int gpuGridSize(int nthreads_per_block);
int gpuMaxNParallel(int idev);
}

#if TINKER_CUDART
#   include <string>
#   include <vector>

namespace tinker {
struct DeviceAttribute
{
   int device;
   std::string name;
   std::string pci_string;

   int cc_major, cc_minor;
   int cc;
   int single_double_ratio;

   std::string compute_mode_string;
   std::string ecc_string;

   size_t free_mem_bytes, total_mem_bytes;

   int max_threads_per_block;
   int max_shared_bytes_per_block;

   int multiprocessor_count;
   int max_threads_per_multiprocessor;
   int max_shared_bytes_per_multiprocessor;
   int max_blocks_per_multiprocessor;
   int cores_per_multiprocessor;

   int clock_rate_kHz; // not memory clock rate
};

std::string gpuCudaRuntimeVersion();
std::string gpuCudaDriverVersion();
std::string gpuThrustVersion();
std::vector<DeviceAttribute>& gpuDeviceAttributes();
}
#endif

#include "mod/gpucard.h"
