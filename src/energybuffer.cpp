#include "tool/energybuffer.h"
#include "ff/atom.h"
#include "tool/darray.h"
#include "tool/rcman.h"
#include <cassert>

namespace tinker {
size_t bufferSize()
{
   size_t bsize = nelem_buffer;
   assert(isPow2(bsize) && "buffer size must be power of 2.");
   return bsize;
}

void bufferAllocate(int flag, EnergyBuffer* pe, VirialBuffer* pv, //
   grad_prec** px, grad_prec** py, grad_prec** pz)
{
   assert(flag & calc::analyz);
   size_t bsize = bufferSize();
   if (flag & calc::energy)
      darray::allocate(bsize, pe);
   if (flag & calc::virial)
      darray::allocate(bsize, pv);
   if (flag & calc::grad)
      darray::allocate(n, px, py, pz);
}

void bufferDeallocate(int flag, EnergyBuffer e, VirialBuffer v, //
   grad_prec* gx, grad_prec* gy, grad_prec* gz)
{
   assert(flag & calc::analyz);
   if (flag & calc::energy)
      darray::deallocate(e);
   if (flag & calc::virial)
      darray::deallocate(v);
   if (flag & calc::grad)
      darray::deallocate(gx, gy, gz);
}

void bufferAllocate(int flag, CountBuffer* pc)
{
   assert(flag & calc::analyz);
   size_t bsize = bufferSize();
   darray::allocate(bsize, pc);
}

void bufferDeallocate(int flag, CountBuffer c)
{
   assert(flag & calc::analyz);
   darray::deallocate(c);
}

int countReduce(const CountBuffer ne)
{
   int c = reduceSum(ne, bufferSize(), g::q0);
   return c;
}

energy_prec energyReduce(const EnergyBuffer e)
{
   auto b = reduceSum(e, bufferSize(), g::q0);
   energy_prec real_out = toFloat<energy_prec>(b);
   return real_out;
}

void virialReduce(virial_prec (&v1)[VirialBufferTraits::N], const VirialBuffer v)
{
   VirialBufferTraits::type b[VirialBufferTraits::N];
   reduceSum2(b, v, bufferSize(), g::q0);
   for (size_t i = 0; i < VirialBufferTraits::N; ++i)
      v1[i] = toFloat<virial_prec>(b[i]);
}

void virialReduce(virial_prec (&v_out)[9], const VirialBuffer v)
{
   virial_prec v1[VirialBufferTraits::N];
   virialReduce(v1, v);
   virialReshape(v_out, v1);
}

void virialReshape(virial_prec (&v_out)[9], const virial_prec (&v1)[VirialBufferTraits::N])
{
   // xx yx zx yy zy zz
   //  0  1  2  3  4  5
   v_out[0] = v1[0]; // xx
   v_out[1] = v1[1]; // xy
   v_out[2] = v1[2]; // xz
   v_out[3] = v1[1]; // yx
   v_out[4] = v1[3]; // yy
   v_out[5] = v1[4]; // yz
   v_out[6] = v1[2]; // zx
   v_out[7] = v1[4]; // zy
   v_out[8] = v1[5]; // zz
}
}
