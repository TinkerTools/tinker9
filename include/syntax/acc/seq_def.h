#pragma once

/**
 * \def SEQ_ROUTINE
 * \ingroup acc_syntax
 * Expands to `_Pragma("acc routine seq")` in OpenACC source files, since we
 * cannot use the `#pragma acc` syntax in a macro.
 */
#define SEQ_ROUTINE _Pragma("acc routine seq")

/**
 * \def SEQ_CUDA
 * \ingroup acc_syntax
 * Is an empty macro in the OpenACC source code.
 */
#define SEQ_CUDA
