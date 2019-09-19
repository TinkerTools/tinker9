#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
static const char* get_SHA1();

void x_info(int argc, char** argv) {
  auto out = stdout;

  print(out, " {}\n\n", "Program Information");

  print(out, "   > This Build is Synchronized with Tinker Commit: {}\n",
        get_SHA1());

#if TINKER_CUDART
  auto loc = "Device with CUDA Runtime Library";
#else
  auto loc = "Host";
#endif
  print(out, "   > Running on {}\n", loc);

  auto rs = sizeof(real);
  print(out, "   > Size of Real: {}\n", sizeof(real));

#if TINKER_DEBUG
  const char* dbg = "ON";
#else
  const char* dbg = "OFF";
#endif
  print(out, "   > Debug: {}\n", dbg);
}

static const char* get_SHA1() {
  return //
      "291a85c1";
  // "6fe8e913";
  // "904bc012";
  // "ddfb803a";
  // "36063480";
}
TINKER_NAMESPACE_END
