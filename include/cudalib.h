#pragma once
#include "error.h"
#include "rc_man.h"
#include <cublas_v2.h>
#include <cuda_runtime.h>


TINKER_NAMESPACE_BEGIN
extern cudaStream_t async_acc;
extern cudaStream_t nonblk;


extern cublasHandle_t h_cublas;
extern cublasHandle_t h_cublas_nonblk;
extern real* pinned_real64;
extern real* dptr_real64;


template <class A, class T1, class T2>
void dotprod(cudaStream_t st, A* ans, const T1* a, const T2* b, int nelem);


void cudalib_data(rc_op);
TINKER_NAMESPACE_END
