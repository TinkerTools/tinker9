#include "tool/cudalib.h"
#include "tool/error.h"

namespace tinker {
void pmeStreamStartRecord(bool usePmeStream)
{
   if (usePmeStream) {
      check_rt(cudaEventRecord(pme_event_start, g::s0));
   }
}

void pmeStreamStartWait(bool usePmeStream)
{
   if (usePmeStream) {
      check_rt(cudaStreamWaitEvent(g::spme, pme_event_start, 0));
   }
}

void pmeStreamFinishRecord(bool usePmeStream)
{
   if (usePmeStream) {
      check_rt(cudaEventRecord(pme_event_finish, g::spme));
   }
}

void pmeStreamFinishWait(bool usePmeStream)
{
   if (usePmeStream) {
      check_rt(cudaStreamWaitEvent(g::s0, pme_event_finish, 0));
   }
}
}
