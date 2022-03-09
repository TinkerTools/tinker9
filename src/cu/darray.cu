#include "tool/error.h"
#include "tool/gpu_card.h"

namespace tinker {
__global__
void device_memory_swap_bytes_async_cu1(char* s1, char* s2, size_t nbytes)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nbytes;
        i += blockDim.x * gridDim.x) {
      char s = s1[i];
      s1[i] = s2[i];
      s2[i] = s;
   }
}

void device_memory_swap_bytes_async_cu(char* s1, char* s2, size_t nbytes,
                                       cudaStream_t st)
{
   int ngrid = get_grid_size(BLOCK_DIM);
   auto ker = device_memory_swap_bytes_async_cu1;
   ker<<<ngrid, BLOCK_DIM, 0, st>>>(s1, s2, nbytes);
   check_rt(cudaGetLastError());
}
}
