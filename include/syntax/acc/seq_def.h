#pragma once


/**
 * \def SEQ_ROUTINE
 * \ingroup macro
 * Expands to `_Pragma("acc routine seq")` in the OpenACC source file.
 * Cannot use the `#pragma` syntax in a macro.
 */
#define SEQ_ROUTINE _Pragma("acc routine seq")


/**
 * \def SEQ_CUDA
 * \ingroup macro
 * Is empty in the OpenACC source code.
 */
#define SEQ_CUDA
