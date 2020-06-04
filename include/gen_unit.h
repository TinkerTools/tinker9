#pragma once
#include "darray.h"
#include <cassert>
#include <memory>


namespace tinker {
enum class GenericUnitVersion
{
   DisableOnDevice,
   EnableOnDevice
};


template <GenericUnitVersion VERS>
struct GenericUnitAlloc;


template <>
struct GenericUnitAlloc<GenericUnitVersion::DisableOnDevice>
{
   static void deallocate(void*) {}
   static void allocate(void**, size_t) {}
   static void copyin(void*, const void*, size_t, DMFlag) {}
};


template <>
struct GenericUnitAlloc<GenericUnitVersion::EnableOnDevice>
{
   static void deallocate(void* p)
   {
      device_memory_deallocate_bytes(p);
   }

   static void allocate(void** pp, size_t nb)
   {
      device_memory_allocate_bytes(pp, nb);
   }

   static void copyin(void* d, const void* s, size_t nb, DMFlag flag)
   {
      device_memory_copyin_bytes(d, s, nb, flag);
   }
};


/**
 * \ingroup rc
 * Resource handle. Analogous to Fortran i/o unit represented by a signed
 * integer.
 */
template <class T,
          GenericUnitVersion VERSION = GenericUnitVersion::DisableOnDevice>
class GenericUnit
{
private:
   static constexpr bool USE_DPTR =
      (VERSION == GenericUnitVersion::DisableOnDevice ? false : true);
   using mem_op = GenericUnitAlloc<VERSION>;
   int unit;


   /**
    * \note
    * The host vector will almost definitely expand its capacity, so if you
    * don't want to implement the move constructors and/or the copy constructors
    * of every possible type T, don't change vector of host pointers to vector
    * of host objects.
    */
   using hostptr_vec = std::vector<std::unique_ptr<T>>;
   using deviceptr_vec =
      std::vector<std::unique_ptr<T, decltype(&mem_op::deallocate)>>;


   static hostptr_vec& hostptrs()
   {
      static hostptr_vec o;
      return o;
   }


   static deviceptr_vec& deviceptrs()
   {
      assert(USE_DPTR);
      static deviceptr_vec o;
      return o;
   }


   const T& obj() const
   {
      assert(0 <= unit && unit < (int)hostptrs().size() &&
             "const T& GenericUnit::obj() const");
      return *hostptrs()[unit];
   }


   T& obj()
   {
      assert(0 <= unit && unit < (int)hostptrs().size() &&
             "T& GenericUnit::obj()");
      return *hostptrs()[unit];
   }


public:
   /// Gets the number of open units.
   static int size()
   {
      if CONSTEXPR (USE_DPTR)
         assert(hostptrs().size() == deviceptrs().size());
      return hostptrs().size();
   }


   /// Releases all of the resources and reset `size()` to 0.
   static void clear()
   {
      // call ~T() on host here
      hostptrs().clear();
      // call deallocate(T*) on device here
      if CONSTEXPR (USE_DPTR)
         deviceptrs().clear();
   }


   /// Resizes the capacity for the objects on host.
   /// \note Cannot be called if device pointers are used.
   template <class DT = T>
   static void resize(int s)
   {
      static_assert(std::is_base_of<T, DT>::value, "");
      assert(!USE_DPTR);
      for (int i = size(); i < s; ++i)
         hostptrs().emplace_back(new DT);
   }


   /// Similar to opening a new Fortran i/o unit.
   /// \return The new unit.
   static GenericUnit open()
   {
      hostptrs().emplace_back(new T);
      if CONSTEXPR (USE_DPTR) {
         T* ptr;
         mem_op::allocate(reinterpret_cast<void**>(&ptr), sizeof(T));
         deviceptrs().emplace_back(ptr, mem_op::deallocate);
      }
      return size() - 1;
   }


   GenericUnit()
      : unit(-1)
   {}


   GenericUnit(int u)
      : unit(u)
   {}


   operator int() const
   {
      return unit;
   }


   /// Whether or not the current unit is open.
   bool valid() const
   {
      return unit >= 0;
   }


   /// Closes the current unit.
   /// \note The resource will not be released until `clear()` is called.
   void close()
   {
      unit = -1;
   }


   /// Gets the (const) reference to the object on host.
   const T& operator*() const
   {
      return obj();
   }
   /// Gets the (const) reference to the object on host.
   T& operator*()
   {
      return obj();
   }


   /// Gets the (const) pointer to the object on host.
   const T* operator->() const
   {
      return &obj();
   }
   /// Gets the (const) pointer to the object on host.
   T* operator->()
   {
      return &obj();
   }


   /// Gets device pointer to the object.
   const T* deviceptr() const
   {
      assert(0 <= unit && unit < deviceptrs().size() &&
             "const T* GenericUnit::deviceptr() const");
      return deviceptrs()[unit].get();
   }
   /// Gets device pointer to the object.
   T* deviceptr()
   {
      assert(0 <= unit && unit < (int)deviceptrs().size() &&
             "T* GenericUnit::deviceptr()");
      return deviceptrs()[unit].get();
   }


   /// Updates the object on device by an object on host.
   /// \param hobj
   /// The reference to the same object on host that can be accessed by the same
   /// unit number.
   void update_deviceptr(const T& hobj, DMFlag flag)
   {
      assert(&hobj == &this->obj());
      mem_op::copyin(this->deviceptr(), &hobj, sizeof(T), flag);
   }
};
}
