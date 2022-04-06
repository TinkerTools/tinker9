#include "tool/ioprint.h"

#include "tinker9.h"

namespace tinker {
void promo()
{
   print(stdout, "%s\n", TINKER9_PROMO_STRING);
}
}
