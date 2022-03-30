#pragma once
#include "tool/rcman.h"
#include <string>
#include <vector>

#if TINKER_CUDART
namespace tinker {
/// \ingroup nvidia
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

/// \ingroup nvidia
std::string gpuCudaRuntimeVersion();
/// \ingroup nvidia
std::string gpuCudaDriverVersion();
/// \ingroup nvidia
std::string gpuThrustVersion();
/// \ingroup nvidia
std::vector<DeviceAttribute>& gpuDeviceAttributes();
}
#endif

namespace tinker {
/// \ingroup platform
void gpuData(RcOp);
/// \ingroup platform
int gpuGridSize(int nthreads_per_block);
/// \ingroup platform
int gpuMaxNParallel(int idev);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup nvidia
/// \brief Number of threads in a warp.
constexpr unsigned WARP_SIZE = 32;
/// \ingroup nvidia
/// \brief Mask for all of the lanes in a warp.
constexpr unsigned ALL_LANES = 0xFFFFFFFF;
/// \ingroup nvidia
/// \brief Default dimension of thread blocks.
constexpr unsigned BLOCK_DIM = /* 64 */ 128 /* 256 */;
/// \ingroup nvidia
constexpr int PME_BLOCKDIM = 64;
/// \ingroup platform
TINKER_EXTERN int ndevice;
/// \ingroup platform
TINKER_EXTERN int idevice;
}
