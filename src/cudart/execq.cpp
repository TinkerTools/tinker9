#include "tool/execq.h"
#include "mod/cudalib.h"
#include "tool/error.h"
#include <cuda_runtime.h>

namespace tinker {
class ExecQ::Impl
{
public:
   cudaStream_t ss;
   cudaEvent_t mdsave_begin_event;
   cudaEvent_t mdsave_end_event;
};

void ExecQ::deallocate()
{
   ptr->ss = nullptr;
   check_rt(cudaEventDestroy(ptr->mdsave_begin_event));
   check_rt(cudaEventDestroy(ptr->mdsave_end_event));
   delete ptr;
}

void ExecQ::allocate()
{
   ptr = new ExecQ::Impl;
   ptr->ss = nullptr;
   check_rt(cudaEventCreateWithFlags(&ptr->mdsave_begin_event, cudaEventDisableTiming));
   check_rt(cudaEventCreateWithFlags(&ptr->mdsave_end_event, cudaEventDisableTiming));
}

void ExecQ::beginCopyout()
{
   check_rt(cudaEventRecord(ptr->mdsave_begin_event, g::s0));
   check_rt(cudaStreamWaitEvent(ptr->ss, ptr->mdsave_begin_event, 0));
}

void ExecQ::endCopyout()
{
   check_rt(cudaEventRecord(ptr->mdsave_end_event, ptr->ss));
   check_rt(cudaStreamWaitEvent(g::s0, ptr->mdsave_end_event, 0));
}
}
