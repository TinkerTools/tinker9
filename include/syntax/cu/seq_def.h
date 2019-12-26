#pragma once


/**
 * \def SEQ_ROUTINE
 * \ingroup macro
 * Expands to `__device__` in the CUDA source file.
 */
#define SEQ_ROUTINE __device__


/**
 * \def SEQ_CUDA
 * \ingroup macro
 * Expands to `__device__` in the CUDA source file.
 * Used in the CUDA kernel templates.
 */
#define SEQ_CUDA __device__
