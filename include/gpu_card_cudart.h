#pragma once
#include "macro.h"
#include <string>
#include <vector>


TINKER_NAMESPACE_BEGIN
struct DeviceAttribute
{
   int device;
   std::string name;

   int cc_major, cc_minor;
   int cc;

   int max_threads_per_block;
   int max_shared_bytes_per_block;

   int multiprocessor_count;
   int max_threads_per_multiprocessor;
   int max_shared_bytes_per_multiprocessor;
   int max_blocks_per_multiprocessor;
};


std::vector<DeviceAttribute>& get_device_attributes();
TINKER_NAMESPACE_END
