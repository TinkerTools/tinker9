#include "intg/baroBerendsen.h"
#include "intg/enum.h"
#include "mdpt.h"

namespace tinker {
BerendsenBarostat::BerendsenBarostat()
   : BasicBarostat(BarostatEnum::Berendsen)
{}

void BerendsenBarostat::control2(time_prec dt)
{
   berendsen_barostat(dt);
}
}
