#include "tool/io_print.h"
#include "version.h"


namespace tinker {
void promo()
{
   print(stdout, "%s\n", TINKER9_PROMO_STRING);
}
}
