#pragma once
#include "mathfunc.h"
#include "mdprec.h"
#include "tool/darray.h"
#include <vector>


namespace tinker {
/**
 * \ingroup mdegv
 * \brief Convert a fixed-point value to floating-point value on host.
 */
template <class T>
inline T to_flt_host(fixed val)
{
   return static_cast<T>(static_cast<long long>(val)) / 0x100000000ull;
}


template <class F, class T>
inline F to_flt_host(T val)
{
   return val;
}


/**
 * \ingroup mdegv
 * \brief
 * The lengths of all of the energy buffers are the same and are implicitly
 * dependent on the number of atoms in the system.
 * \note Must be a power of 2.
 */
size_t buffer_size();


/**
 * \ingroup mdegv
 * \brief
 * Poor man's buffer array for energy, virial tensor, etc. accumulation.
 * In fact, the buffer itself is not fancy. It is just a raw pointer to a
 * pre-allocated array on device. This class stores some properties that
 * cannot easily be fetched from the raw pointers.
 *
 * This table shows a few possible definitions of the buffers (which may or
 * may not be the definitions actually used).
 *
 * |                | `T[N]`         | `type (*)[value]` | `type[value]`       |
 * |----------------|----------------|-------------------|---------------------|
 * | **Properties** | **Increments** | **Pointer Types** | **Buffer Elements** |
 * | Counts         | int            | int*              | int                 |
 * | Virial         | double[9]      | double (*)[16]    | double[16]          |
 * | Virial (symm.) | float[6]       | fixed (*)[8]      | fixed[8]            |
 * | Energy         | float          | fixed*            | fixed               |
 *
 * \tparam T      Type of energy increments.
 * \tparam Nincr  Number of elements per increment; 1 in general;
 *                6 for symmetric virial tensor; 9 for general virial tensor.
 *
 * The class (at least) consists of following members:
 *    - type[value]: Type of underlying buffer elements; `value` must be
 *    a power of 2 and must not be not less than Nincr; by default, `T` and
 *    `type` are the same.
 *    - N: Must be equal to Nicr.
 */
template <class T, size_t Nincr>
struct buffer_traits
{
   static constexpr size_t N = Nincr;
   static constexpr size_t value = pow2_ge(Nincr);
   using type = T;
};


/**
 * \ingroup mdegv
 * \brief Special rules for `float` increments where fixed-point buffer is used.
 */
template <size_t Nincr>
struct buffer_traits<float, Nincr>
{
   static constexpr size_t N = Nincr;
   static constexpr size_t value = pow2_ge(Nincr);
   using type = fixed;
};


using count_buffer_traits = buffer_traits<int, 1>;
using energy_buffer_traits = buffer_traits<e_prec, 1>;
using virial_buffer_traits = buffer_traits<v_prec, 6>;


//====================================================================//


using count_buffer =
   pointer<count_buffer_traits::type, count_buffer_traits::value>;
using energy_buffer =
   pointer<energy_buffer_traits::type, energy_buffer_traits::value>;
using virial_buffer =
   pointer<virial_buffer_traits::type, virial_buffer_traits::value>;


void buffer_allocate(int flag, energy_buffer*, virial_buffer*, grad_prec**,
                     grad_prec**, grad_prec**);
void buffer_deallocate(int flag, energy_buffer, virial_buffer, grad_prec*,
                       grad_prec*, grad_prec*);
void buffer_allocate(int, count_buffer*);
void buffer_deallocate(int, count_buffer);


/**
 * \ingroup mdegv
 * \brief
 * Get the number of the non-bonded interactions, energy, or virial from the
 * corresponding buffer. These operations unnecessarily finish within `O(1)`
 * time complexity.
 */
int count_reduce(const count_buffer b);
/**
 * \ingroup mdegv
 * \see count_reduce
 */
energy_prec energy_reduce(const energy_buffer b);
/**
 * \ingroup mdegv
 * \see count_reduce
 */
void virial_reduce(virial_prec (&)[virial_buffer_traits::N],
                   const virial_buffer b);
/**
 * \ingroup mdegv
 * \see count_reduce
 */
void virial_reduce(virial_prec (&)[9], const virial_buffer b);
}
