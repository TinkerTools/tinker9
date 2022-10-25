// ck.py Version 3.1.0
__global__
static void mdPos_cu1(int n, time_prec dt, pos_prec* restrict qx, pos_prec* restrict qy, pos_prec* restrict qz,
   const vel_prec* restrict vlx, const vel_prec* restrict vly, const vel_prec* restrict vlz)
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      qx[i] += dt * vlx[i];
      qy[i] += dt * vly[i];
      qz[i] += dt * vlz[i];
   }
}
