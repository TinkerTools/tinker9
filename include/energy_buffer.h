#pragma once
#include "darray.h"
#include "mathfunc.h"
#include "mdprec.h"
#include <vector>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup md_egv
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
 * \ingroup md_egv
 * \brief
 * The lengths of all of the energy buffers are the same and are implicitly
 * dependent on the number of atoms in the system.
 * \note
 * Must a power of 2.
 */
size_t buffer_size();


/**
 * \ingroup md_egv
 * \brief
 * Poor man's buffer array for energy, virial tensor, etc. accumulation.
 * In fact, the buffer itself is not fancy. It is just a raw pointer to a
 * pre-allocated array on device. This class defines some properties that
 * cannot easily be stored by or fetched from the raw pointers.
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
 * \ingroup md_egv
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


/**
 * \ingroup md_egv
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
void buffer_allocate(energy_buffer*, grad_prec**, grad_prec**, grad_prec**,
                     virial_buffer*);
/**
 * \ingroup md_egv
 * \see buffer_allocate
 */
void buffer_deallocate(energy_buffer, grad_prec*, grad_prec*, grad_prec*,
                       virial_buffer);
/**
 * \ingroup md_egv
 * \see buffer_allocate
 */
void buffer_allocate(count_buffer*);
/**
 * \ingroup md_egv
 * \see buffer_allocate
 */
void buffer_deallocate(count_buffer);


/**
 * \ingroup md_egv
 * \brief
 * Get the number of the non-bonded interactions, energy, or virial from the
 * corresponding buffer. These operations unnecessarily finish within `O(1)`
 * time complexity.
 */
int count_reduce(const count_buffer b);
/**
 * \ingroup md_egv
 * \see count_reduce
 */
energy_prec energy_reduce(const energy_buffer b);
/**
 * \ingroup md_egv
 * \see count_reduce
 */
void virial_reduce(virial_prec (&)[virial_buffer_traits::N],
                   const virial_buffer b);
/**
 * \ingroup md_egv
 * \see count_reduce
 */
void virial_reduce(virial_prec (&)[9], const virial_buffer b);


/**
 * \ingroup md_egv
 * \brief Bookkeeping list for the count buffers.
 */
extern std::vector<count_buffer> count_buffers;
/**
 * \ingroup md_egv
 * \brief Bookkeeping list for the energy buffers.
 */
extern std::vector<energy_buffer> energy_buffers;
/**
 * \ingroup md_egv
 * \brief Bookkeeping list for the virial buffers.
 */
extern std::vector<virial_buffer> virial_buffers;
/**
 * \ingroup md_egv
 */
extern std::vector<grad_prec*> x_grads, y_grads, z_grads;
TINKER_NAMESPACE_END
