#include "tool/gpu_card.h"
#include "algorithm"
#include "md.h"
#include "tinker_rt.h"
#include "tool/error.h"
#include "tool/exec.h"
#include <cuda_runtime.h>
#include <thrust/version.h>


namespace tinker {
std::string get_cuda_runtime_version_string()
{
   int ver, major, minor;
   check_rt(cudaRuntimeGetVersion(&ver));
   // ver = 1000*major + 10*minor
   major = ver / 1000;
   minor = (ver - major * 1000) / 10;
   return format("%d.%d", major, minor);
}


std::string get_cuda_driver_version_string()
{
   int ver, major, minor;
   check_rt(cudaDriverGetVersion(&ver));
   // ver = 1000*major + 10*minor
   major = ver / 1000;
   minor = (ver - major * 1000) / 10;
   return format("%d.%d", major, minor);
}


std::string get_thrust_version_string()
{
   return format("%d.%d.%d patch %d", THRUST_MAJOR_VERSION,
                 THRUST_MINOR_VERSION, THRUST_SUBMINOR_VERSION,
                 THRUST_PATCH_NUMBER);
}


std::vector<DeviceAttribute>& get_device_attributes()
{
   static std::vector<DeviceAttribute> a;
   return a;
}


static std::string get_nvidia_smi()
{
   std::string smi = "nvidia-smi";
   int val1 = std::system("which nvidia-smi > /dev/null");
   if (val1 != 0) {
      val1 = std::system("which nvidia-smi.exe > /dev/null");
      smi = "nvidia-smi.exe";
      if (val1 != 0) {
         TINKER_THROW("nvidia-smi is not found.");
      }
   }
   return smi;
}


static void get_device_attribute(DeviceAttribute& a, int device = 0)
{
   cudaDeviceProp prop;
   check_rt(cudaGetDeviceProperties(&prop, device));


   a.device = device;
   a.name = prop.name;


   char pciBusID[16];
   check_rt(cudaDeviceGetPCIBusId(pciBusID, 13, device));
   a.pci_string = pciBusID;


   a.cc_major = prop.major;
   a.cc_minor = prop.minor;
   a.cc = a.cc_major * 10 + a.cc_minor;
   a.single_double_ratio = prop.singleToDoublePrecisionPerfRatio;


   if (prop.computeMode == cudaComputeModeExclusive)
      a.compute_mode_string = "exclusive thread";
   else if (prop.computeMode == cudaComputeModeProhibited)
      a.compute_mode_string = "prohibited";
   else if (prop.computeMode == cudaComputeModeExclusiveProcess)
      a.compute_mode_string = "exclusive process";
   else
      a.compute_mode_string = "default";
   if (prop.ECCEnabled)
      a.ecc_string = "on";
   else
      a.ecc_string = "off";


   check_rt(cudaSetDevice(device));
   check_rt(cudaMemGetInfo(&a.free_mem_bytes, &a.total_mem_bytes));


   a.max_threads_per_block = prop.maxThreadsPerBlock;
   a.max_shared_bytes_per_block = prop.sharedMemPerBlock;

   a.multiprocessor_count = prop.multiProcessorCount;
   a.max_threads_per_multiprocessor = prop.maxThreadsPerMultiProcessor;
   a.max_shared_bytes_per_multiprocessor = prop.sharedMemPerMultiprocessor;


   // nvidia cuda-c-programming-guide compute-capabilities
   bool found_cc = true;


   // Maximum number of resident blocks per multiprocessor
   if (a.cc > 75)
      found_cc = false;
   else if (a.cc >= 75)
      a.max_blocks_per_multiprocessor = 16;
   else if (a.cc >= 50)
      a.max_blocks_per_multiprocessor = 32;
   else if (a.cc >= 30)
      a.max_blocks_per_multiprocessor = 16;
   else
      found_cc = false;


   // Number of CUDA cores per multiprocessor, not tabulated in
   // cuda-c-programming-guide;
   // documented in "Compute Capability - architecture"
   // 7.0 7.2 7.5: 64
   // 6.1 6.2: 128
   // 6.0: 64
   // 5.0 5.2: 128
   // 3.0 3.5 3.7: 192
   if (a.cc > 75)
      found_cc = false;
   else if (a.cc >= 70)
      a.cores_per_multiprocessor = 64;
   else if (a.cc >= 61)
      a.cores_per_multiprocessor = 128;
   else if (a.cc >= 60)
      a.cores_per_multiprocessor = 64;
   else if (a.cc >= 50)
      a.cores_per_multiprocessor = 128;
   else if (a.cc >= 30)
      a.cores_per_multiprocessor = 192;
   else
      found_cc = false;


   a.clock_rate_kHz = prop.clockRate;


   check_rt(cudaDeviceReset());


   if (!found_cc) {
      TINKER_THROW(
         format("The code base should be updated for compute capability %d; "
                "Please refer to the Nvidia Cuda-C Programming Guide",
                a.cc));
   }
}


static int recommend_device(int ndev)
{
   int usp = -1; // user-specified cuda device; -1 for not set
   const char* usp_str = nullptr;


   // if do not use xyz file, there is no key file
   if (usp < 0 && (rc_flag & calc::xyz)) {
      usp_str = "CUDA-DEVICE keyword";
      get_kv("CUDA-DEVICE", usp, -1);
   }
   // check environment variable "CUDA_DEVICE"
   if (usp < 0) {
      usp_str = "CUDA_DEVICE environment variable";
      if (const char* str = std::getenv("CUDA_DEVICE"))
         usp = std::stoi(str);
   }
   // check environment variable "cuda_device"
   if (usp < 0) {
      usp_str = "cuda_device environment variable";
      if (const char* str = std::getenv("cuda_device"))
         usp = std::stoi(str);
   }


   std::vector<int> gpercent, prcd; // precedence
   std::vector<double> gflops;
   for (int i = 0; i < ndev; ++i) {
      const auto& a = get_device_attributes()[i];
      std::string smi = get_nvidia_smi();
      std::string cmd = format("%s --query-gpu=utilization.gpu "
                               "--format=csv,noheader,nounits -i %s",
                               smi, a.pci_string);
      std::string percent = exec(cmd);
      prcd.push_back(i);
      gpercent.push_back(std::stoi(percent));
      double gf = a.clock_rate_kHz;
      gf *= a.cores_per_multiprocessor * a.multiprocessor_count;
      gflops.push_back(gf);
   }


   auto strictly_prefer = [&](int idev, int jdev) {
      // If idev is preferred return true.
      // If idev and jdev are considered equivalent return false.
      // If jdev is preferred return false.
      const int SIGNIFICANT_DIFFERENCE = 8;


      int igpuutil = gpercent[idev];
      int jgpuutil = gpercent[jdev];
      // choose the idle device
      // if the difference is significant, stop here
      if (igpuutil + SIGNIFICANT_DIFFERENCE < jgpuutil)
         return true;
      else if (jgpuutil + SIGNIFICANT_DIFFERENCE < igpuutil)
         return false;


      double igflp = gflops[idev];
      double jgflp = gflops[jdev];
      // choose the faster device
      // if the difference is significant, stop here
      if (igflp > jgflp * (1.0 + SIGNIFICANT_DIFFERENCE / 100.0))
         return true;
      else if (jgflp > igflp * (1.0 + SIGNIFICANT_DIFFERENCE / 100.0))
         return false;


      // If idev is preferred, the program would never reach this line.
      return false;
   };
   std::stable_sort(prcd.begin(), prcd.end(), strictly_prefer);


   int idev;
   if (usp < 0) {
      usp_str = "GPU utilization";
      idev = prcd[0];
   } else if (usp == prcd[0]) {
      idev = usp;
   } else {
      print(stdout,
            "\n CUDA-DEVICE Warning,"
            " Program recommended Device %d but Device %d was set from %s\n",
            prcd[0], usp, usp_str);
      idev = usp;
   }
   print(stdout, "\n GPU Device :  Setting Device ID to %d from %s\n", idev,
         usp_str);
   return idev;
}


static unsigned int cuda_device_flags = 0;
void gpu_card_data(rc_op op)
{
   if (op & rc_dealloc) {
      /*
      // should not reset these variables for the unit tests, if there are
      // multiple GPUs available
      ndevice = 0;

      get_device_attributes().clear();

      idevice = -1;
      */
   }

   if (op & rc_init) {
      if (cuda_device_flags) {
         return;
      }

      always_check_rt(cudaGetDeviceCount(&ndevice));

      auto& all = get_device_attributes();
      all.resize(ndevice);
      for (int i = 0; i < ndevice; ++i)
         get_device_attribute(all[i], i);

      idevice = recommend_device(ndevice);

      check_rt(cudaSetDevice(idevice));

      if (cuda_device_flags == 0) {
         cuda_device_flags = cudaDeviceMapHost;
#if 1
         cuda_device_flags |= cudaDeviceScheduleBlockingSync;
#elif 0
         // Using this flag may reduce the latency
         // for cudaStreamSynchronize() calls.
         cuda_device_flags |= cudaDeviceScheduleSpin;
#endif
         always_check_rt(cudaSetDeviceFlags(cuda_device_flags));
      }
   }
}


int get_grid_size(int nthreads_per_block)
{
   const auto& a = get_device_attributes()[idevice];

   nthreads_per_block = std::min(nthreads_per_block, a.max_threads_per_block);
   int max_nblocks_per_MP =
      std::min((a.max_threads_per_multiprocessor + nthreads_per_block - 1) /
                  nthreads_per_block,
               a.max_blocks_per_multiprocessor);

   return a.multiprocessor_count * max_nblocks_per_MP;
}


int gpu_max_nparallel(int idev)
{
   const auto& a = get_device_attributes().at(idev);
   int n_sm = a.multiprocessor_count;
   int n_thread = a.max_threads_per_multiprocessor;
   return n_sm * n_thread;
}
}
