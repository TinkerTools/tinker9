#include <CoreFoundation/CoreFoundation.h>
#include <IOKit/graphics/IOGraphicsLib.h>

int tinkerGpuUtilizationInt32_macos(int iDevice)
{
   CFMutableDictionaryRef dict = IOServiceMatching(kIOAcceleratorClassName);

   io_iterator_t it;
   if (IOServiceGetMatchingServices(kIOMasterPortDefault, dict, &it) == kIOReturnSuccess) {

      int i = 0;
      CFMutableDictionaryRef serviceDictionary;
      io_registry_entry_t entry;
      while ((entry = IOIteratorNext(it))) {
         if (IORegistryEntryCreateCFProperties(
                entry, &serviceDictionary, kCFAllocatorDefault, kNilOptions) != kIOReturnSuccess) {
            IOObjectRelease(entry);
            return -1; // failed to create service dictionary
         }

         ssize_t gpuCoreUse = 0;
         CFMutableDictionaryRef perfProperties = (CFMutableDictionaryRef)CFDictionaryGetValue(
            serviceDictionary, CFSTR("PerformanceStatistics"));
         if (perfProperties) {
            const void* gpuCoreUtilization =
               CFDictionaryGetValue(perfProperties, CFSTR("GPU Core Utilization"));
            if (gpuCoreUtilization)
               CFNumberGetValue((CFNumberRef)gpuCoreUtilization, kCFNumberSInt64Type, &gpuCoreUse);
         }

         CFRelease(serviceDictionary);
         IOObjectRelease(entry);

         if (i == iDevice) {
            IOObjectRelease(it);
            return gpuCoreUse / 10000000;
         }
         ++i;
      } // end iterate over the device found

      IOObjectRelease(it);
   }

   return -1;
}
