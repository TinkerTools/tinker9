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
 * rc_man foo_random_name{foo_data, rc};
 * rc_man bar_random_name{bar_data, rc};
 * @endcode
 *
 * To deallocate resource in the same order of allocation, use unnamed objects.
 * @code
 * rc_man{foo_data, rc};
 * rc_man{bar_data, rc};
 * @endcode
 */
class rc_man {
public:
  typedef enum {
    rc_dealloc = 0x001, /// deallocate device memory
    rc_alloc = 0x002,   /// allocate device memory
    rc_copyin = 0x004,  /// update device data from host memory, or directly
                        /// initialize device data
    rc_evolve = 0x010,

    dealloc = rc_dealloc,
    alloc = rc_alloc,
    init = rc_copyin,
    evolve = rc_evolve
  } rc_t;

private:
  void (*f_)(rc_t);
  rc_t rc_;
  bool dealloc_() const { return rc_ & dealloc; }

public:
  rc_man(void (*f)(rc_t), rc_t rc) : f_(f), rc_(rc) {
    if (!dealloc_()) {
      f_(rc_);
    }
  }

  ~rc_man() {
    if (dealloc_()) {
      f_(rc_);
    }
  }
};
typedef rc_man::rc_t rc_t;
constexpr rc_t rc_dealloc = rc_t::dealloc;
constexpr rc_t rc_alloc = rc_t::alloc;
constexpr rc_t rc_copyin = rc_t::init;
constexpr rc_t rc_evolve = rc_t::evolve;
typedef rc_t rc_op;

void host_data(rc_op);
void device_data(rc_op);
TINKER_NAMESPACE_END

#include "util_cudart.h"
#include "util_hostonly.h"
#include "util_io.h"

#define m_tinker_get_3rd_arg_(arg1, arg2, arg3, ...) arg3
#define m_tinker_get_4th_arg_(arg1, arg2, arg3, arg4, ...) arg4

#define m_tinker_always_check_cudart_1_(cucall)                                \
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

#define m_tinker_always_check_cudart_2_(cucall, optmsg)                        \
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

#define m_tinker_always_check_cudart_3_(cucall, res_t, cu_0)                   \
  do {                                                                         \
    res_t cures_ = cucall;                                                     \
    if (cures_ != cu_0) {                                                      \
      print_backtrace();                                                       \
      std::string m_ = format(" errno {} of type {} at {}:{}", cures_,         \
                              TINKER_STR(res_t), __FILE__, __LINE__);          \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define m_tinker_always_check_cudart_(...)                                     \
  m_tinker_get_4th_arg_(__VA_ARGS__, m_tinker_always_check_cudart_3_,          \
                        m_tinker_always_check_cudart_2_,                       \
                        m_tinker_always_check_cudart_1_)

#if TINKER_DEBUG || defined(TINKER_ALWAYS_CHECK_CUDART)
#  define check_cudart(...)                                                    \
    m_tinker_always_check_cudart_(__VA_ARGS__)(__VA_ARGS__)
#else
#  define check_cudart(cucall, ...) cucall
#endif

#endif
