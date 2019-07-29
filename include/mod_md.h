#ifndef TINKER_MOD_MD_H_
#define TINKER_MOD_MD_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
typedef enum {
  thermo_berendsen,
  thermo_bussi,
  thermo_andersen,
  thermo_nose_hoover_chain,
  thermo_null
} thermostat_t;
TINKER_EXTERN thermostat_t thermostat;

typedef enum {
  baro_berendsen,
  baro_bussi,
  baro_nose_hoover_chain,
  baro_montecarlo,
  baro_null
} barostat_t;
TINKER_EXTERN barostat_t barostat;
}
TINKER_NAMESPACE_END

#endif
