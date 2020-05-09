#pragma once
#include "darray.h"
#include "mathfunc.h"
#include "mdprec.h"
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


/**
 * \ingroup mdegv
 * \brief
 * Allocate or deallocate device memory for the buffers.
 *
 * ### Count Buffers
 * Only be allocated only if #calc::analyz flag is used.
 *
 *
 * ### Energy
 *
 *
 * ### Virial
 *
 *
 * ### Gradients
 *
 *
 * Using bond and VDW terms as examples:
 *
 * |                 | if analyz  | if .not. analyz |
 * |-----------------|------------|-----------------|
 * | total potential | #esum      | #esum           |
 * | total gx        | #gx        | #gx             |
 * | total virial    | #vir       | #vir            |
 * |                 |            |                 |
 * | ebond buffer    | #eb        | #eng_buf        |
 * | ebond           | #energy_eb | #esum           |
 * | ebond gx        | #debx      | #gx             |
 * | vir_eb buffer   | #vir_buf   | #vir_buf        |
 * | vir_eb          | #vir       | #vir            |
 * |                 |            |                 |
 * | evdw buffer     | #ev        | #ev             |
 * | evdw            | #energy_ev | #energy_ev      |
 * | evdw gx         | #devx      | #devx           |
 * | vir_ev buffer   | #vir_buf   | #vir_buf        |
 * | vir_ev          | #vir       | #vir            |
 *
 * \note
 * There is a global list to bookkeep all of the allocated buffers of its own
 * kind, so that all of the buffers can be iterated.
 * \note
 * These functions cannot be used on #eng_buf and #vir_buf.
 */
void buffer_allocate(int rcflag, energy_buffer* ebuf, grad_prec** gx,
                     grad_prec** gy, grad_prec** gz, virial_buffer* vbuf,
                     energy_prec* eout);
void buffer_deallocate(int, energy_buffer, grad_prec*, grad_prec*, grad_prec*,
                       virial_buffer);
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


/**
 * \ingroup mdegv
 */
void set_energy_reduce_dst(energy_buffer, energy_prec*);
/**
 * \ingroup mdegv
 */
energy_prec* get_energy_reduce_dst(energy_buffer);
/**
 * \ingroup mdegv
 */
void clear_energy_reduce_dst();


/**
 * \ingroup mdegv
 * \brief Bookkeeping list for the count buffers.
 */
extern std::vector<count_buffer> count_buffers;
/**
 * \ingroup mdegv
 * \brief Bookkeeping list for the energy buffers.
 */
extern std::vector<energy_buffer> energy_buffers;
/**
 * \ingroup mdegv
 * \brief Bookkeeping list for the virial buffers.
 */
extern std::vector<virial_buffer> virial_buffers;
/**
 * \ingroup mdegv
 * \brief Bookkeeping list for the gradients.
 */
extern std::vector<grad_prec*> x_grads;
/**
 * \ingroup mdegv
 * \copydoc x_grads
 */
extern std::vector<grad_prec*> y_grads;
/**
 * \ingroup mdegv
 * \copydoc x_grads
 */
extern std::vector<grad_prec*> z_grads;
}
