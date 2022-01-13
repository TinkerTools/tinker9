#include "intg/baroBerendsen.h"
#include "mdpt.h"

namespace tinker {
void BerendsenBarostat::control2(time_prec dt)
{
   berendsen_barostat(dt);
}
}
