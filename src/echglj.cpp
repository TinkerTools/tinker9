#include "echglj.h"
#include "nblist.h"
#include "potent.h"


namespace tinker {
void echglj_data(rc_op op)
{
   if (!(clist_version() & NBL_SPATIAL))
      return;


   if (!use_potent(charge_term) || !use_potent(vdw_term))
      return;
}


void echglj(int vers) {}
}
