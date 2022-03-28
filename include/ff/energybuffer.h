#pragma once
#include "math/pow2.h"
#include "precision.h"
#include "tool/darray.h"

namespace tinker {
/// \brief Convert a fixed-point value to floating-point value on host.
template <class T>
inline T toFloat_host(fixed val)
{
   return static_cast<T>(static_cast<long long>(val)) / 0x100000000ull;
}

template <class F, class T>
inline F toFloat_host(T val)
{
   return val;
}

/// \brief The lengths of all of the energy buffers are the same and are implicitly
/// dependent on the number of atoms in the system.
/// \note Must be a power of 2.
size_t bufferSize();

/// \brief
/// Poor man's buffer array for energy, virial tensor, etc. accumulation.
/// In fact, the buffer itself is not fancy. It is just a raw pointer to a
/// pre-allocated array on device. This class stores some properties that
/// cannot easily be fetched from the raw pointers.
///
/// This table shows a few possible definitions of the buffers (which may or
/// may not be the definitions actually used).
///
/// |                | `T[N]`         | `type (*)[value]` | `type[value]`       |
/// |----------------|----------------|-------------------|---------------------|
/// | **Properties** | **Increments** | **Pointer Types** | **Buffer Elements** |
/// | Counts         | int            | int*              | int                 |
/// | Virial         | double[9]      | double (*)[16]    | double[16]          |
/// | Virial (symm.) | float[6]       | fixed (*)[8]      | fixed[8]            |
/// | Energy         | float          | fixed*            | fixed               |
///
/// \tparam T      Type of energy increments.
/// \tparam Nincr  Number of elements per increment; 1 in general;
///                6 for symmetric virial tensor; 9 for general virial tensor.
///
/// The class (at least) consists of following members:
///    - type[value]: Type of underlying buffer elements; `value` must be
///    a power of 2 and must not be not less than Nincr; by default, `T` and
///    `type` are the same.
///    - N: Must be equal to Nicr.
template <class T, size_t Nincr>
struct BufferTraits
{
   static constexpr size_t N = Nincr;
   static constexpr size_t value = pow2Ge(Nincr);
   using type = T;
};

/// \brief Special rules for `float` increments where fixed-point buffer is used.
template <size_t Nincr>
struct BufferTraits<float, Nincr>
{
   static constexpr size_t N = Nincr;
   static constexpr size_t value = pow2Ge(Nincr);
   using type = fixed;
};

using CountBufferTraits = BufferTraits<int, 1>;
using EnergyBufferTraits = BufferTraits<e_prec, 1>;
using VirialBufferTraits = BufferTraits<v_prec, 6>;

using CountBuffer = pointer<CountBufferTraits::type, CountBufferTraits::value>;
using EnergyBuffer = pointer<EnergyBufferTraits::type, EnergyBufferTraits::value>;
using VirialBuffer = pointer<VirialBufferTraits::type, VirialBufferTraits::value>;

void bufferAllocate(int flag, EnergyBuffer*, VirialBuffer*, grad_prec**, grad_prec**, grad_prec**);
void bufferDeallocate(int flag, EnergyBuffer, VirialBuffer, grad_prec*, grad_prec*, grad_prec*);
void bufferAllocate(int, CountBuffer*);
void bufferDeallocate(int, CountBuffer);

/// \brief
/// Get the number of the non-bonded interactions, energy, or virial from the
/// corresponding buffer. These operations unnecessarily finish within `O(1)`
/// time complexity.
int countReduce(const CountBuffer b);
/// \see countReduce
energy_prec energyReduce(const EnergyBuffer b);
void virialReshape(virial_prec (&output)[9], const virial_prec (&input)[VirialBufferTraits::N]);
/// \see countReduce
void virialReduce(virial_prec (&)[VirialBufferTraits::N], const VirialBuffer b);
/// \see countReduce
void virialReduce(virial_prec (&)[9], const VirialBuffer b);
}
