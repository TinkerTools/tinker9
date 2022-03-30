#pragma once

/// \def SEQ_ROUTINE
/// \ingroup acc_syntax
/// \brief Expands to `_Pragma("acc routine seq")` in OpenACC source files.
/// `#pragma acc` cannot be used in macro.
#define SEQ_ROUTINE _Pragma("acc routine seq")

/// \def SEQ_CUDA
/// \ingroup acc_syntax
/// \brief An empty macro in the OpenACC source code.
#define SEQ_CUDA
