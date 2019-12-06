#include "gpu_card.h"
#include "error.h"
#include "mathfunc.h"
#include <cuda_runtime.h>
#include <limits>
#include <map>
#include <thrust/version.h>


TINKER_NAMESPACE_BEGIN
std::string get_cuda_runtime_version_string()
{
   int ver, major, minor;
   check_rt(cudaRuntimeGetVersion(&ver));
   // ver = 1000*major + 10*minor
   major = ver / 1000;
   minor = (ver - major * 1000) / 10;
   return format("{}.{}", major, minor);
}


std::string get_cuda_driver_version_string()
{
   int ver, major, minor;
   check_rt(cudaDriverGetVersion(&ver));
   // ver = 1000*major + 10*minor
   major = ver / 1000;
   minor = (ver - major * 1000) / 10;
   return format("{}.{}", major, minor);
}


std::string get_thrust_version_string()
{
   return format("{}.{}.{} patch {}", THRUST_MAJOR_VERSION,
                 THRUST_MINOR_VERSION, THRUST_SUBMINOR_VERSION,
                 THRUST_PATCH_NUMBER);
}


std::vector<DeviceAttribute>& get_device_attributes()
{
   static std::vector<DeviceAttribute> a;
   return a;
}


static void get_device_attribute(DeviceAttribute& a, int device = 0)
{
   cudaDeviceProp prop;
   check_rt(cudaGetDeviceProperties(&prop, device));

   a.device = device;
   a.name = prop.name;

   a.pci_string = format("{:02x}:{:02x}.{}", prop.pciBusID, prop.pciDeviceID,
                         prop.pciDomainID);


   a.cc_major = prop.major;
   a.cc_minor = prop.minor;
   a.cc = a.cc_major * 10 + a.cc_minor;
   a.single_double_ratio = prop.singleToDoublePrecisionPerfRatio;


   if (prop.computeMode == cudaComputeModeExclusive)
      a.compute_mode_string = "Exclusive Thread";
   else if (prop.computeMode == cudaComputeModeProhibited)
      a.compute_mode_string = "Prohibited";
   else if (prop.computeMode == cudaComputeModeExclusiveProcess)
      a.compute_mode_string = "Exclusive Process";
   else
      a.compute_mode_string = "Default";
   if (prop.ECCEnabled)
      a.ecc_string = "ON";
   else
      a.ecc_string = "OFF";


   check_rt(cudaMemGetInfo(&a.free_mem_bytes, &a.total_mem_bytes));


   a.max_threads_per_block = prop.maxThreadsPerBlock;
   a.max_shared_bytes_per_block = prop.sharedMemPerBlock;

   a.multiprocessor_count = prop.multiProcessorCount;
   a.max_threads_per_multiprocessor = prop.maxThreadsPerMultiProcessor;
   a.max_shared_bytes_per_multiprocessor = prop.sharedMemPerMultiprocessor;

   // nvida cuda-c-programming-guide compute-capabilities

   int value;

   value = 0;
   switch (a.cc) {
   case 30:
   case 32:
   case 35:
   case 37:
      a.max_blocks_per_multiprocessor = 16;
      break;
   case 50:
   case 52:
   case 53:
   case 60:
   case 61:
   case 62:
   case 70:
      a.max_blocks_per_multiprocessor = 32;
      break;
   case 75:
      a.max_blocks_per_multiprocessor = 16;
      break;
   default:
      a.max_blocks_per_multiprocessor = value;
      break;
   }
   if (a.max_blocks_per_multiprocessor <= 0) {
      TINKER_THROW(
         format("Cannot get maximum number of resident blocks per "
                "multiprocessor for the undocumented compute capability {}; "
                "Please refer to the Nvida Cuda-C Programming Guide",
                a.cc));
   }
}


static unsigned int cuda_device_flags = 0;
void gpu_card_data(rc_op op)
{
   if (op & rc_dealloc) {
      ndevice = 0;

      get_device_attributes().clear();

      idevice = -1;
   }

   if (op & rc_init) {
      if (cuda_device_flags == 0) {
         cuda_device_flags = cudaDeviceMapHost;
         cuda_device_flags |= cudaDeviceScheduleBlockingSync;
         check_rt(cudaSetDeviceFlags(cuda_device_flags));
      }

      check_rt(cudaGetDeviceCount(&ndevice));

      auto& all = get_device_attributes();
      all.resize(ndevice);
      for (int i = 0; i < ndevice; ++i)
         get_device_attribute(all[i], i);

      idevice = 0;
      check_rt(cudaSetDevice(idevice));
   }
}


int get_grid_size(int nthreads_per_block)
{
   const auto& a = get_device_attributes()[idevice];

   nthreads_per_block = std::min(nthreads_per_block, a.max_threads_per_block);
   int max_nblocks_per_MP =
      std::min(a.max_threads_per_multiprocessor / nthreads_per_block,
               a.max_blocks_per_multiprocessor);

   return a.multiprocessor_count * max_nblocks_per_MP;
}


static int get_block_size1(int shared_bytes_per_thread)
{
   const auto& a = get_device_attributes()[idevice];
   const int bsize_min = WARP_SIZE;

   // determine max block size

   int bsize_max = a.max_threads_per_block;
   int upper_bound = std::numeric_limits<int>::max();
   if (shared_bytes_per_thread > 0) {
      upper_bound = a.max_shared_bytes_per_block / shared_bytes_per_thread;
      bsize_max = std::min(bsize_max, upper_bound);
   }

   int bsize = 0;
   int occ = 0;
   const int max_warps_mproc = a.max_threads_per_multiprocessor / WARP_SIZE;
   for (int b = bsize_min; b <= bsize_max; b *= 2) {
      int max_blocks_per_mproc = a.max_blocks_per_multiprocessor;
      int upper_bound = std::numeric_limits<int>::max();
      if (shared_bytes_per_thread > 0) {
         upper_bound = a.max_shared_bytes_per_multiprocessor /
            (b * shared_bytes_per_thread);
         max_blocks_per_mproc = std::min(max_blocks_per_mproc, upper_bound);
      }

      // calculate occupancy with this block size

      int warps_per_block = (b + WARP_SIZE - 1) / WARP_SIZE;
      int actual_nblocks_per_mproc = max_warps_mproc / warps_per_block;
      actual_nblocks_per_mproc =
         std::min(actual_nblocks_per_mproc, max_blocks_per_mproc);

      int occup = warps_per_block * actual_nblocks_per_mproc;
      if (occup > occ) {
         bsize = b;
         occ = occup;
      }
   }

   if (bsize <= 0) {
      TINKER_THROW(
         format("The value of Shared Bytes per Thread ({}) is too large",
                shared_bytes_per_thread));
   }

   assert(ct::is_pow2(bsize));
   return bsize;
}


int get_block_size(int shared_bytes_per_thread)
{
   struct Config
   {
      int device;
      int shared_bytes_per_thread;

      bool operator<(const Config& c) const
      {
         if (device == c.device) {
            return shared_bytes_per_thread < c.shared_bytes_per_thread;
         } else {
            return device < c.device;
         }
      }
   };
   static std::map<Config, int> config_to_bsize;

   shared_bytes_per_thread = std::max(0, shared_bytes_per_thread);

   Config c;
   c.device = idevice;
   c.shared_bytes_per_thread = shared_bytes_per_thread;
   auto it = config_to_bsize.find(c);
   if (it == config_to_bsize.end()) {
      int s = get_block_size1(shared_bytes_per_thread);
      config_to_bsize[c] = s;
      return s;
   } else {
      return it->second;
   }
}
TINKER_NAMESPACE_END
