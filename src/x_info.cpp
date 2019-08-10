#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
static const char* get_SHA1();

void x_info(int argc, char** argv) {
  auto out = stdout;

  print(out, " {}\n\n", "Program Information");

  print(out, "   > This Build is Synchronized with Tinker Commit: {}\n",
        get_SHA1());

#if defined(TINKER_CUDART)
  auto loc = "Device with CUDA Runtime Library";
#else
  auto loc = "Host";
#endif
  print(out, "   > Running on {}\n", loc);

  auto rs = sizeof(real);
  print(out, "   > Size of Real: {}\n", sizeof(real));
}

static const char* get_SHA1() {
  return //
      "291a85c1435feddc835e80bfa340497b67cc1393";
  // "6fe8e913fe4da3d46849d10248ad2a4872b4da93";
  // "904bc0125aa0ab548866bc3effed34df1ec1b4d6";
  // "ddfb803a9a35237dd624a9eb9e74ddde4ea9062f";
  // "36063480f115bd0455d4d6291f9666ab5f0b739d";
}
TINKER_NAMESPACE_END
