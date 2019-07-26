#include "gpu/gpudata.h"
#include "gpu/gpu.h"
#include "util_rc_man.h"

extern "C" {
void tinker_gpu_data_create() {
  using namespace TINKER_NAMESPACE;
  using namespace gpu;
  const rc_t rc = static_cast<rc_t>(rc_alloc | rc_copyin);
  host_data(rc);
  gpu_data(rc);
}

void tinker_gpu_data_create_() { tinker_gpu_data_create_(); }

void tinker_gpu_data_destroy() {
  using namespace TINKER_NAMESPACE;
  using namespace gpu;
  const rc_t rc = rc_dealloc;
  gpu_data(rc);
  host_data(rc);
}

void tinker_gpu_data_destroy_() { tinker_gpu_data_destroy(); }
}

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern void n_data(rc_t);

extern void xyz_data(rc_t);
extern void vel_data(rc_t);
extern void mass_data(rc_t);

extern void potential_data(rc_t);

extern void box_data(rc_t);
extern void couple_data(rc_t);
extern void nblist_data(rc_t);

extern void md_data(rc_t);

void gpu_data(rc_t rc) {
  rc_man n42_{n_data, rc};

  rc_man xyz42_{xyz_data, rc};
  rc_man vel42_{vel_data, rc};
  rc_man mass42_{mass_data, rc};

  rc_man pd42_{potential_data, rc};

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
  rc_man box42_{box_data, rc};
  rc_man cpl42_{couple_data, rc};
  rc_man nbl42_{nblist_data, rc};

  rc_man md42_{md_data, rc};
}
}
TINKER_NAMESPACE_END
