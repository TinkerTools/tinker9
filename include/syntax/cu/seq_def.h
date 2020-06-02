#pragma once


/**
 * \def SEQ_ROUTINE
 * \ingroup cuda_syntax
 * Expands to `__device__` in CUDA source files.
 */
#define SEQ_ROUTINE __device__


/**
 * \def SEQ_CUDA
 * \ingroup cuda_syntax
 * Expands to `__device__` in CUDA source files.
 * Used in CUDA kernel templates.
 */
#define SEQ_CUDA __device__
