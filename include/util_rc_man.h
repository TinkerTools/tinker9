#ifndef TINKER_UTIL_RC_MAN_H_
#define TINKER_UTIL_RC_MAN_H_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN
void fortran_runtime_initialize(int, char**);
void fortran_runtime_finish();
void tinker_gpu_runtime_initialize();
void tinker_gpu_runtime_finish();

/**
 * Resource Management
 *
 * To deallocate resource in reverse order of allocation, use named objects.
 * @code
 * rc_man foo_random_name{foo_data, op};
 * rc_man bar_random_name{bar_data, op};
 * @endcode
 *
 * To deallocate resource in the same order of allocation, use unnamed objects.
 * @code
 * rc_man{foo_data, op};
 * rc_man{bar_data, op};
 * @endcode
 */
class ResourceManagement {
public:
  typedef enum {
    dealloc = 0x001, ///< deallocate device memory
    alloc = 0x002,   ///< allocate device memory
    init = 0x004,    ///< update device data from host memory, or directly
                     ///< initialize device data

    evolve = 0x010
  } op_t;

private:
  void (*f_)(op_t);
  op_t op_;
  bool will_dealloc_() const { return op_ & dealloc; }

public:
  ResourceManagement(void (*f)(op_t), op_t op)
      : f_(f)
      , op_(op) {
    if (!will_dealloc_()) {
      f_(op_);
    }
  }

  ~ResourceManagement() {
    if (will_dealloc_()) {
      f_(op_);
    }
  }
};
typedef ResourceManagement rc_man;
typedef rc_man::op_t rc_op;
constexpr rc_op rc_dealloc = rc_man::dealloc;
constexpr rc_op rc_alloc = rc_man::alloc;
constexpr rc_op rc_init = rc_man::init;

void host_data(rc_op);
void device_data(rc_op);
TINKER_NAMESPACE_END

#include "util_cudart.h"
#include "util_hostonly.h"
#include "util_io.h"

#define TINKER_GET_3RD_ARG_(arg1, arg2, arg3, ...) arg3
#define TINKER_GET_4TH_ARG_(arg1, arg2, arg3, arg4, ...) arg4

#define TINKER_ALWAYS_CHECK_CUDART_1_(cucall)                                  \
  do {                                                                         \
    cudaError_t cures_ = cucall;                                               \
    if (cures_ != cudaSuccess) {                                               \
      print_backtrace();                                                       \
      const char* msg = cudaGetErrorString(cures_);                            \
      std::string m_ =                                                         \
          format(" {} (errno {}) at {}:{}", msg, cures_, __FILE__, __LINE__);  \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define TINKER_ALWAYS_CHECK_CUDART_2_(cucall, optmsg)                          \
  do {                                                                         \
    cudaError_t cures_ = cucall;                                               \
    if (cures_ != cudaSuccess) {                                               \
      print_backtrace();                                                       \
      const char* msg = cudaGetErrorName(cures_);                              \
      std::string m_ = format(" {} {} (errno {}) at {}:{}", optmsg, msg,       \
                              cures_, __FILE__, __LINE__);                     \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define TINKER_ALWAYS_CHECK_CUDART_3_(cucall, res_t, cu_0)                     \
  do {                                                                         \
    res_t cures_ = cucall;                                                     \
    if (cures_ != cu_0) {                                                      \
      print_backtrace();                                                       \
      std::string m_ = format(" errno {} of type {} at {}:{}", cures_,         \
                              TINKER_STR(res_t), __FILE__, __LINE__);          \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define TINKER_ALWAYS_CHECK_CUDART_(...)                                       \
  TINKER_GET_4TH_ARG_(__VA_ARGS__, TINKER_ALWAYS_CHECK_CUDART_3_,              \
                      TINKER_ALWAYS_CHECK_CUDART_2_,                           \
                      TINKER_ALWAYS_CHECK_CUDART_1_)

#if TINKER_DEBUG || defined(TINKER_ALWAYS_CHECK_CUDART)
#  define check_rt(...) TINKER_ALWAYS_CHECK_CUDART_(__VA_ARGS__)(__VA_ARGS__)
#else
#  define check_rt(cucall, ...) cucall
#endif

#endif
