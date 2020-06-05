#include "tinker_rt.h"
#include "tool/compilers.h"
#include "tool/gpu_card.h"
#include "version.h"


namespace tinker {
namespace {
const char* get_SHA1()
{
   return         //
      "350df099"; // Tue Apr  7 21:01:20 2020 -0500
   // "11e84c69"; // Tue Nov 12 14:56:12 2019 -0600
   // "291a85c1"; // Fri Jul 19 16:21:27 2019 +0200
   // "6fe8e913"; // Sun Apr 21 13:34:28 2019 -0500
   // "904bc012";
   // "ddfb803a";
   // "36063480";
}


std::string get_version_string()
{
   std::string r = format("%d.%d.%d", TINKER_GPU_VERSION_MAJOR,
                          TINKER_GPU_VERSION_MINOR, TINKER_GPU_VERSION_PATCH);
#ifdef TINKER_GPU_GIT_SHORT_HASH
   r += format(" GIT %s", TINKER_STR(TINKER_GPU_GIT_SHORT_HASH));
#endif
   return r;
}
}

void x_info(int, char**)
{
   auto out = stdout;
   auto fmt = "    %-36s %s\n";
   auto fmz = "    %-36s %zu\n";


   platform_data(rc_init);
   gpu_card_data(rc_init);


   print(out, "\n\n Program Information\n\n");
   print(out, fmt, "Version:", get_version_string());
   print(out, fmt, "Synchronized with Tinker commit:", get_SHA1());
   print(out, fmt, "C++ compiler:", cxx_compiler_name());
   print(out, fmz, "Size of real (bytes):", sizeof(real));
   print(out, fmz, "Size of mixed (bytes):", sizeof(mixed));
   print(out, fmt, "Using deterministic force:",
         TINKER_DETERMINISTIC_FORCE ? "true" : "false");
#if TINKER_DEBUG
   const char* dbg = "on";
#else
   const char* dbg = "off";
#endif
   print(out, fmt, "Debug mode:", dbg);


#if TINKER_HOST
   print(out, fmt, "Platform:", "CPU");
#endif


#if TINKER_CUDART
   auto fmd = "    %-36s %d\n";
   auto fm1 = "    %-36s\n";
   auto fm2 = "       %-33s %s\n";
   auto f2d = "       %-33s %d\n";


   print(out, fmt, "Platform:", "CUDA and OpenACC");
   if (pltfm_config & CU_PLTFM)
      print(out, fmt, "Primary GPU package:", "CUDA");
   else if (pltfm_config & ACC_PLTFM)
      print(out, fmt, "Primary GPU package:", "OpenACC");
   print(out, fmt,
         "Latest CUDA supported by driver:", get_cuda_driver_version_string());
   print(out, fmt, "CUDA runtime version:", get_cuda_runtime_version_string());
   print(out, fmt, "Thrust version:", get_thrust_version_string());
   print(out, fmt, "CUDA compiler:", cuda_compiler_name());
   print(out, fmt, "OpenACC compiler:", acc_compiler_name());
   if (ndevice > 0) {
      print(out, fmd, "GPU detected:", ndevice);
      const auto& attribs = get_device_attributes();
      for (const auto& a : attribs) {
         print(out, fm1, format("GPU %d:", a.device));
         print(out, fm2, "PCI:", a.pci_string);
         print(out, fm2, "Name:", a.name);
         print(out, fm2, "Maximum compute capability:",
               format("%d.%d", a.cc_major, a.cc_minor));
         print(out, f2d, "Single double perf. ratio:", a.single_double_ratio);
         print(out, fm2, "Compute mode:", a.compute_mode_string);
         print(out, fm2, "Error-correcting code (ECC):", a.ecc_string);
         print(out, f2d, "Clock rate (kHz):", a.clock_rate_kHz);
         print(out, f2d, "Number of Multiprocessors:", a.multiprocessor_count);
         print(out, f2d, "Number of CUDA cores:",
               a.cores_per_multiprocessor * a.multiprocessor_count);
         const double B_to_GB = 1024. * 1024. * 1024.;
         print(out, fm2, "Used/Total GPU memory:",
               format("%.2f % / %.2f GB",
                      100.0 - 100.0 * a.free_mem_bytes / a.total_mem_bytes,
                      a.total_mem_bytes / B_to_GB));
      }
   }
#endif


   gpu_card_data(rc_dealloc);
   platform_data(rc_dealloc);
}
}
