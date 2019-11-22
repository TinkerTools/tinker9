#pragma once
#include "dev_array.h"
#include <vector>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup ebuf
 * \brief
 * The lengths of all of the energy buffers are the same and are implicitly
 * dependent on the number of atoms in the system. Depending on the template
 * parameter `N`, the underlying buffer is either equivalent to
 * `type[buffer_size()]` if `N=1`, or `type[value][buffer_size()]` if `N>1`.
 * \see buffer_traits::n buffer_traits::value buffer_traits::type
 */
size_t buffer_size();


/**
 * \ingroup ebuf
 * \brief
 * A variable (an energy, a virial, etc.) of type `T[N]` or `T` (if `N`
 * equals 1) can be stored in a buffer to get better performance and higher
 * accuracy.
 *
 * Value `N` is also used as `n` inside the class. The type of the underlying
 * buffer elements is `type[value]`, which may or may not be the same as `T[N]`.
 * If the types are different, this class provides a cast function to convert
 * `type` back to `T`, otherwise, this cast function does nothing. `value` is a
 * power of 2 and is not less than `n`.
 *
 * This class provides access to the properties described above.
 */
template <class T, size_t N>
struct buffer_traits
{
   static constexpr size_t n = N;
   static constexpr size_t value = ct::pow2_ge(N);
   typedef T type;
   static T cast(type val)
   {
      return val;
   }
};


template <size_t N>
struct buffer_traits<float, N>
{
   static constexpr size_t n = N;
   static constexpr size_t value = ct::pow2_ge(N);
   typedef unsigned long long type;
   static float cast(type val)
   {
      return static_cast<float>(static_cast<long long>(val)) / 0x100000000ull;
   }
};


/// \ingroup ebuf
/// \see buffer_traits device_pointer
using count_buffer_traits = buffer_traits<int, 1>;
/// \ingroup ebuf
/// \see buffer_traits device_pointer
using energy_buffer_traits = buffer_traits<real, 1>;
/// \ingroup ebuf
/// \see buffer_traits device_pointer
using virial_buffer_traits = buffer_traits<real, 6>;
/// \ingroup ebuf
/// \see buffer_traits device_pointer
using count_buffer =
   device_pointer<count_buffer_traits::type, count_buffer_traits::value>;
/// \ingroup ebuf
/// \see buffer_traits device_pointer
using energy_buffer =
   device_pointer<energy_buffer_traits::type, energy_buffer_traits::value>;
/// \ingroup ebuf
/// \see buffer_traits device_pointer
using virial_buffer =
   device_pointer<virial_buffer_traits::type, virial_buffer_traits::value>;


/**
 * \ingroup ebuf
 * \brief
 * Allocate or deallocate device memory for the buffers.
 * \note
 * There is a global list to bookkeep all of the allocated buffers of its own
 * kind, so that all of the buffers can be iterated.
 * \note
 * These functions cannot be used on `esum_buf` and `vir_buf`.
 * \note
 * If `calc::analyz` is not used, the allocate functions only point the energy
 * buffer and virial buffer to `esum_buf` and `vir_buf`, respectively, and point
 * the count buffer to `nullptr`. If `calc::analyz` is used, `esm_buf` and
 * `vir_buf` are all set to `nullptr`.
 * \see count_buffers energy_buffers virial_buffers
 * \see esum_buf vir_buf
 * \see calc::analyz
 */
void buffer_allocate(energy_buffer*, virial_buffer*);
/// \ingroup ebuf
/// \see buffer_allocate
void buffer_deallocate(energy_buffer, virial_buffer);
/// \ingroup ebuf
/// \see buffer_allocate
void buffer_allocate(count_buffer*, energy_buffer*, virial_buffer*);
/// \ingroup ebuf
/// \see buffer_allocate
void buffer_deallocate(count_buffer, energy_buffer, virial_buffer);


/**
 * \ingroup ebuf
 * \brief
 * Get the number of the non-bonded interactions, energy, or virial from the
 * corresponding buffer. These operations unnecessarily finish within `O(1)`
 * time complexity.
 */
int get_count(const count_buffer b);
/// \ingroup ebuf
/// \see get_count
real get_energy(const energy_buffer b);
/// \ingroup ebuf
/// \see get_count
void get_virial(real (&)[virial_buffer_traits::n], const virial_buffer b);
/// \ingroup ebuf
/// \see get_count
void get_virial(real (&)[9], const virial_buffer b);


/// \ingroup ebuf
/// \brief Bookkeeping list for the count buffers.
extern std::vector<count_buffer> count_buffers;
/// \ingroup ebuf
/// \brief Bookkeeping list for the energy buffers.
extern std::vector<energy_buffer> energy_buffers;
/// \ingroup ebuf
/// \brief Bookkeeping list for the virial buffers.
extern std::vector<virial_buffer> virial_buffers;
TINKER_NAMESPACE_END
