#ifndef TINKER_ENERGY_BUFFER_H_
#define TINKER_ENERGY_BUFFER_H_

#include "dev_array.h"
#include "gen_unit.h"
#include "mathfunc.h"

TINKER_NAMESPACE_BEGIN
namespace ebuf_detail {
/// @brief
/// convert fixed point scalar @c val to a floating point number
template <class T>
struct FixedToFloating {
  static T exec(FixedPointType val) {
    assert(std::is_floating_point<T>::value);
    return static_cast<T>(static_cast<long long>(val)) / fixed_point;
  }
};

/// @brief
/// @c S<T>::Type gives the underlying type of the buffer; by default using the
/// same type as @c T; use partial specialization for special cases
/// @{
template <class T>
struct S {
  typedef T Type;
};

template <>
struct S<float> {
  typedef FixedPointType Type;
};

template <>
struct S<double> {
  typedef double Type;
  // typedef FixedPointType Type;
};
/// @}

/// @brief
/// sum device array @c dptr[@c nelem][@c NStore] and save to
/// host array @c host_ans[@c Answer]
template <class Answer, class Store, int NAnswer, int NStore>
struct Sum {
  static void exec(Answer* __restrict__ host_ans,
                   const Store* __restrict__ dptr, int nelem) {
    static_assert(NAnswer >= 1, "");
    static_assert(NStore >= NAnswer, "");

    Store val[NAnswer];
    if_constexpr(NAnswer == 1) { val[0] = reduce_sum(dptr, nelem); }
    else {
      reduce_sum2(val, NAnswer, dptr, nelem, NStore);
    }

    if_constexpr(std::is_same<Store, FixedPointType>::value &&
                 !std::is_same<Store, Answer>::value) {
      if_constexpr(NAnswer == 1) {
        *host_ans = FixedToFloating<Answer>::exec(val[0]);
      }
      else {
        for (int i = 0; i < NAnswer; ++i)
          host_ans[i] = FixedToFloating<Answer>::exec(val[i]);
      }
    }
    else if_constexpr(std::is_same<Store, Answer>::value) {
      if_constexpr(NAnswer == 1) { *host_ans = val[0]; }
      else {
        for (int i = 0; i < NAnswer; ++i)
          host_ans[i] = val[i];
      }
    }
    else {
      assert(false);
    }
  }
};
}

//====================================================================//

/// @brief
/// each element of the buffer on device is of type @c Store[@c NStore];
/// can be reduced to type @c Answer[@c NAnswer] on host;
/// the length of the underlying device array must be a power of 2
template <class Answer, int NAnswer, int NStore>
class GenericBuffer {
private:
  static_assert(NStore >= NAnswer, "");
  typedef typename ebuf_detail::S<Answer>::Type Store;
  static constexpr size_t RS = sizeof(Store) * NStore;

  Store* buf;
  int cap;

  void grow_if_must(int new_size) {
    if (new_size <= cap)
      return;

    int old_cap = cap;
    int new_cap = new_size;

    Store* new_buf;
    device_array::allocate(RS * new_cap, &new_buf);
    device_array::copy(RS * old_cap, new_buf, buf);
    device_array::deallocate(buf);

    buf = new_buf;
    cap = new_cap;
  }

public:
  static int calc_size(int nelem, int max_MB = 4) {
    size_t max_bytes = max_MB * 1024 * 1024ull;
    size_t new_size = max_bytes / sizeof(Answer);

    assert(is_pow2(new_size) && "new_size must be power of 2");
    return new_size;
  }

public:
  typedef Store* PointerType;
  const Store* buffer() const { return buf; }
  Store* buffer() { return buf; }
  int size() const { return cap; }
  void zero() { device_array::zero(RS * cap, buf); }
  void sum(Answer* host_ans) {
    ebuf_detail::Sum<Answer, Store, NAnswer, NStore>::exec(host_ans, buf, cap);
  }

  void alloc(int nelem, int max_MB = 4) {
    int new_size = calc_size(nelem, max_MB);
    grow_if_must(new_size);
  }

  GenericBuffer()
      : buf(nullptr)
      , cap(0) {}

  ~GenericBuffer() {
    cap = 0;
    device_array::deallocate(buf);
    buf = nullptr;
  }
};

//====================================================================//

typedef GenericBuffer<int, 1, 1> CountBuffer;
typedef GenericBuffer<real, 1, 1> EnergyBuffer;
typedef GenericBuffer<real, 9, 16> VirialBuffer;

typedef GenericUnit<CountBuffer, GenericUnitVersion::DisableOnDevice> Count;
typedef GenericUnit<EnergyBuffer, GenericUnitVersion::DisableOnDevice> Energy;
typedef GenericUnit<VirialBuffer, GenericUnitVersion::DisableOnDevice> Virial;

int get_count(Count);
double get_energy(Energy);
void get_virial(double*, Virial);

class BondedEnergy {
protected:
  int bufsize_;
  Energy e_;
  Virial vir_;

public:
  int buffer_size() const;
  Energy e();
  Virial vir();

  void dealloc();
  void alloc(int bsize);
};

class NonbondedEnergy : public BondedEnergy {
protected:
  Count ne_;

public:
  Count ne();

  void dealloc();
  void alloc(int bsize);
};
TINKER_NAMESPACE_END

#endif
