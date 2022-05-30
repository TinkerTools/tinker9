#pragma once
#include "tool/rcman.h"
#include <string>
#include <vector>

#if TINKER_CUDART
namespace tinker {
/// \addtogroup nvidia
/// \{

/// Device attributes
struct DeviceAttribute
{
   int device;             ///< Device number.
   std::string name;       ///< Device name.
   std::string pci_string; ///< PCI Bus ID.

   int cc_major;            ///< Major compute capability.
   int cc_minor;            ///< Minor compute capability.
   int cc;                  ///< Compute capability multiplied by 10.
   int single_double_ratio; ///< Single to double precision performance ratio.

   std::string
      compute_mode_string; ///< `exclusive thread`, `prohibited`, `exclusive process`, or `default`.
   std::string ecc_string; ///< `on`, or `off`.

   size_t free_mem_bytes;  ///< Free device memory in bytes.
   size_t total_mem_bytes; ///< Total device memory in bytes.

   int max_threads_per_block;      ///< Maximum threads per block.
   int max_shared_bytes_per_block; ///< Maximum shared memory per block in bytes.

   int multiprocessor_count;                ///< Multiprocessor count.
   int max_threads_per_multiprocessor;      ///< Maximum threads per multiprocessor.
   int max_shared_bytes_per_multiprocessor; ///< Maximum shared memory per mutiprocessor in bytes.
   int max_blocks_per_multiprocessor;       ///< Maximum thread blocks per multiporcessor.
   int cores_per_multiprocessor;            ///< Number of cores per multiprocessor.

   int clock_rate_kHz; ///< Clock frequency in kHz, not memory clock rate.
};

/// \return  CUDA runtime version.
std::string gpuCudaRuntimeVersion();
/// \return  Max CUDA runtime version supported by the driver.
std::string gpuCudaDriverVersion();
/// \return  Version of the Thrust Library.
std::string gpuThrustVersion();
/// \return  Attributes of all visible devices.
std::vector<DeviceAttribute>& gpuDeviceAttributes();

/// \}
}
#endif

namespace tinker {
/// \addtogroup platform
/// \{

/// Sets up the GPU card.
void gpuData(RcOp);

/// \param nthreadsPerBlock  Dimension of the thread blocks.
/// \return                  Recommended dimension of the thread grids.
int gpuGridSize(int nthreadsPerBlock);

/// \param idev  Device number.
/// \return      Maximum number of parallel threads on the device.
int gpuMaxNParallel(int idev);

///\}
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \addtogroup nvidia
/// \{
constexpr unsigned WARP_SIZE = 32;         ///< Number of threads in a warp.
constexpr unsigned ALL_LANES = 0xFFFFFFFF; ///< Mask for all of the lanes in a warp.
constexpr unsigned BLOCK_DIM = 128;        ///< Default dimension of thread blocks.
constexpr int PME_BLOCKDIM = 64;           ///< Dimension of the thread blocks for some PME kernels.
/// \}

/// \addtogroup platform
/// \{
TINKER_EXTERN int ndevice; ///< Total number of visible devices.
TINKER_EXTERN int idevice; ///< Device ID number in use.
/// \}
}
