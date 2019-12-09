#include "e_mpole.h"


TINKER_NAMESPACE_BEGIN
void empole_coulomb_cu(int vers)
{
   extern void empole_coulomb_acc(int);
   empole_coulomb_acc(vers);
}
TINKER_NAMESPACE_END
