#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
void x_info(int argc, char** argv) {
  auto out = stdout;

  print(out, " {}\n\n", "Program Information");

#ifdef TINKER_HOST
  auto loc = "Host";
#else
  auto loc = "Device";
#endif
  print(out, "   > Running on {}\n", loc);

  auto rs = sizeof(real);
  print(out, "   > Size of Real: {}\n", sizeof(real));
}
TINKER_NAMESPACE_END
