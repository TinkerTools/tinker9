#pragma once
#include "macro.h"
#include <string>
#include <vector>


TINKER_NAMESPACE_BEGIN
std::string get_cuda_runtime_version_string();
std::string get_cuda_driver_version_string();
std::string get_thrust_version_string();


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
};


std::vector<DeviceAttribute>& get_device_attributes();
TINKER_NAMESPACE_END
