#include "gpu_card.h"
#include "launch.h"
#include "md.h"
#include "syntax/cu/reduce.h"
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>


namespace tinker {
// velocity to ekin[3][3] (actually ekin[6])
template <unsigned int B>
__global__
void velocity_to_ekin_cu(energy_prec* out, const vel_prec* restrict vx,
                         const vel_prec* restrict vy,
                         const vel_prec* restrict vz,
                         const mass_prec* restrict mass, int n,
                         energy_prec ekcal_inv)
{
   constexpr int HN = 6;
   __shared__ energy_prec sd[HN][B];
   unsigned int t = threadIdx.x;
   #pragma unroll
   for (int j = 0; j < HN; ++j)
      sd[j][t] = 0;
   for (int i = t + blockIdx.x * B; i < n; i += B * gridDim.x) {
      energy_prec term = 0.5f * mass[i] * ekcal_inv;
      sd[0][t] += term * vx[i] * vx[i]; // exx
      sd[1][t] += term * vy[i] * vy[i]; // eyy
      sd[2][t] += term * vz[i] * vz[i]; // ezz
      sd[3][t] += term * vx[i] * vy[i]; // exy
      sd[4][t] += term * vy[i] * vz[i]; // eyz
      sd[5][t] += term * vz[i] * vx[i]; // ezx
   }
   __syncthreads();


   using Op = OpPlus<energy_prec>;
   Op op;
   // clang-format off
   if (B >= 512) { if (t < 256) { _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 256]); } __syncthreads(); }
   if (B >= 256) { if (t < 128) { _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 128]); } __syncthreads(); }
   if (B >= 128) { if (t < 64 ) { _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 64 ]); } __syncthreads(); }
   if (t < 32  ) warp_reduce2<energy_prec, HN, B, Op>(sd, t, op);
   // clang-format on
   if (t == 0)
      #pragma unroll
      for (int j = 0; j < HN; ++j)
         out[blockIdx.x * HN + j] = sd[j][0];
}


void kinetic_cu(T_prec& temp)
{
   cudaStream_t st = nonblk;
   constexpr int HN = 6;
   energy_prec* dptr = reinterpret_cast<energy_prec*>(dptr_buf);
   energy_prec(*dptr6)[HN] = reinterpret_cast<energy_prec(*)[HN]>(dptr_buf);
   energy_prec* hptr = reinterpret_cast<energy_prec*>(pinned_buf);
   int grid_siz1 = get_grid_size(BLOCK_DIM);
   grid_siz1 = grid_siz1 / HN;
   int grid_siz2 = (n + BLOCK_DIM - 1) / BLOCK_DIM;
   int grid_size = std::min(grid_siz1, grid_siz2);
   const energy_prec ekcal_inv = 1.0 / units::ekcal;
   velocity_to_ekin_cu<BLOCK_DIM>
      <<<grid_size, BLOCK_DIM, 0, st>>>(dptr, vx, vy, vz, mass, n, ekcal_inv);
   reduce2<energy_prec, BLOCK_DIM, HN, HN, OpPlus<energy_prec>>
      <<<1, BLOCK_DIM, 0, st>>>(dptr6, dptr6, grid_size);
   check_rt(cudaMemcpyAsync(hptr, dptr, HN * sizeof(energy_prec),
                            cudaMemcpyDeviceToHost, st));
   check_rt(cudaStreamSynchronize(st));
   energy_prec exx = hptr[0];
   energy_prec eyy = hptr[1];
   energy_prec ezz = hptr[2];
   energy_prec exy = hptr[3];
   energy_prec eyz = hptr[4];
   energy_prec ezx = hptr[5];


   ekin[0][0] = exx;
   ekin[0][1] = exy;
   ekin[0][2] = ezx;
   ekin[1][0] = exy;
   ekin[1][1] = eyy;
   ekin[1][2] = eyz;
   ekin[2][0] = ezx;
   ekin[2][1] = eyz;
   ekin[2][2] = ezz;
   eksum = exx + eyy + ezz;
   temp = 2 * eksum / (mdstuf::nfree * units::gasconst);
}
}
