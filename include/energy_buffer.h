#ifndef TINKER_ENERGY_BUFFER_H_
#define TINKER_ENERGY_BUFFER_H_

#include "dev_array.h"
#include "gen_unit.h"
#include "mathfunc.h"

TINKER_NAMESPACE_BEGIN
/**
 * \brief
 * A buffer class to store the "energy" components, where the total energy can
 * be obtained from the reduction operation, although unnecessarily within
 * `O(1)` time.
 *
 * \tparam TA,NA
 * \c TA and \c NA are referred to as \c A and \c N, respectively.
 * The type of the total "energy" is \c A if \c N is 1,
 * or `A[N]` if \c N is greater than 1.
 * Type \c A also dictates the underlying storing type \c S on device.
 * For examples, number of interactions: `int` or `size_t`;
 * single precision energy `float`;
 * double precision virial `double[6]` or `double[9]`, etc.
 *
 * \note
 * The length of the underlying device array must be a power of 2.
 */
template <class TA, size_t NA>
class GenericBuffer
{
public:
   typedef TA A;
   static constexpr size_t N = NA;

   //====================================================================//

private:
   /// \brief
   /// `Type` gives the underlying type of the device array.
   template <class DUMMY, class T>
   struct S0
   {
      typedef T Type;
   };

   template <class DUMMY>
   struct S0<DUMMY, float>
   {
      typedef fixed Type;
   };

   template <class DUMMY>
   struct S0<DUMMY, double>
   {
      typedef double Type;
      // typedef fixed Type;
   };

public:
   typedef typename S0<void, A>::Type S;
   static constexpr size_t NS = pow2_ge (N);

private:
   static_assert (NS >= N, "");
   static_assert (is_pow2 (NS), "");
   static constexpr size_t RS = sizeof (S) * NS;

   //====================================================================//

   /// \brief
   /// Convert a fixed point scalar \c val to a floating point number.
   template <class T>
   static T fixed_to_floating (fixed val)
   {
      assert (std::is_floating_point<T>::value);
      return static_cast<T> (static_cast<long long> (val)) / fixed_point;
   }

   //====================================================================//

   /// \brief
   /// Sum device array `dptr[nelem][NS]` and save to host array
   /// `host_ans[A]`.
   template <class Store, size_t NS>
   static void sum (A* __restrict__ host_ans, const Store* __restrict__ dptr,
                    size_t nelem)
   {
      static_assert (N >= 1, "");
      static_assert (NS >= N, "");

      Store val[N];
      if_constexpr (N == 1)
      {
         val[0] = parallel::reduce_sum (dptr, nelem);
      }
      else
      {
         parallel::reduce_sum2 (val, N, dptr, nelem, NS);
      }

      if_constexpr (std::is_same<Store, fixed>::value &&
                    !std::is_same<Store, A>::value)
      {
         if_constexpr (N == 1)
         {
            *host_ans = fixed_to_floating<A> (val[0]);
         }
         else
         {
            for (size_t i = 0; i < N; ++i)
               host_ans[i] = fixed_to_floating<A> (val[i]);
         }
      }
      else if_constexpr (std::is_same<Store, A>::value)
      {
         if_constexpr (N == 1)
         {
            *host_ans = val[0];
         }
         else
         {
            for (size_t i = 0; i < N; ++i)
               host_ans[i] = val[i];
         }
      }
      else
      {
         assert (false);
      }
   }

   //====================================================================//

   S* buf_;
   size_t cap_;

   void grow_if_must (size_t new_size)
   {
      if (new_size <= cap_)
         return;

      size_t old_cap = cap_;
      size_t new_cap = new_size;

      S* new_buf;
      device_array::allocate (RS * new_cap, &new_buf);
      device_array::copy (RS * old_cap, new_buf, buf_);
      device_array::deallocate (buf_);

      buf_ = new_buf;
      cap_ = new_cap;
   }

public:
   static size_t calc_size (size_t nelem)
   {
      size_t max_bytes = 4 * 1024 * 1024ull; // 4 MB for word size
      max_bytes *= (sizeof (A) / sizeof (int));
      if (nelem <= 16384)
         max_bytes /= 2; // 2 MB
      if (nelem <= 8192)
         max_bytes /= 2; // 1 MB
      if (nelem <= 4096)
         max_bytes /= 2; // 512 KB
      if (nelem <= 2048)
         max_bytes /= 2; // 256 KB
      if (nelem <= 1024)
         max_bytes /= 2; // 128 KB
      if (nelem <= 512)
         max_bytes /= 2; // 64 KB 16384 words

      size_t new_size = max_bytes / sizeof (A);

      assert (is_pow2 (new_size) && "new_size must be power of 2");
      return new_size;
   }

public:
   const S* buffer () const
   {
      return buf_;
   }
   S* buffer ()
   {
      return buf_;
   }
   size_t size () const
   {
      return cap_;
   }
   void zero ()
   {
      device_array::zero (RS * cap_, buf_);
   }
   void sum (A* host_ans)
   {
      sum<S, NS> (host_ans, buf_, cap_);
   }

   void alloc (size_t nelem)
   {
      size_t new_size = calc_size (nelem);
      grow_if_must (new_size);
   }

   GenericBuffer ()
      : buf_ (nullptr)
      , cap_ (0)
   {}

   ~GenericBuffer ()
   {
      cap_ = 0;
      device_array::deallocate (buf_);
      buf_ = nullptr;
   }
};

//====================================================================//

typedef GenericBuffer<int, 1> CountBuffer;
typedef GenericBuffer<real, 1> EnergyBuffer;
typedef GenericBuffer<real, 6> VirialBuffer;

typedef GenericUnit<CountBuffer, GenericUnitVersion::DisableOnDevice> Count;
typedef GenericUnit<EnergyBuffer, GenericUnitVersion::DisableOnDevice> Energy;
typedef GenericUnit<VirialBuffer, GenericUnitVersion::DisableOnDevice> Virial;

int get_count (Count);
double get_energy (Energy);
void get_virial (double*, Virial);

class BondedEnergy
{
protected:
   size_t bufsize_;
   Energy e_;
   Virial vir_;

public:
   size_t buffer_size () const;
   Energy e ();
   Virial vir ();

   void dealloc ();
   void alloc (size_t bsize);
};

class NonbondedEnergy : public BondedEnergy
{
protected:
   Count ne_;

public:
   Count ne ();

   void dealloc ();
   void alloc (size_t bsize);
};
TINKER_NAMESPACE_END

#endif
