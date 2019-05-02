#ifndef TINKER_GPU_DATA_H_
#define TINKER_GPU_DATA_H_

#include "defines.h"

extern "C" {
void tinker_gpu_data_create();
void tinker_gpu_data_destroy();
void tinker_gpu_gradient();
}

#endif
