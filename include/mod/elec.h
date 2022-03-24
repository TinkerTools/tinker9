#pragma once
#include "ff/energybuffer.h"
#include "ff/pme.h"
#include "precision.h"

namespace tinker {
TINKER_EXTERN real electric;
TINKER_EXTERN real dielec;
}

// pme
namespace tinker {
TINKER_EXTERN PMEUnit epme_unit;  // electrostatic
TINKER_EXTERN PMEUnit ppme_unit;  // polarization
TINKER_EXTERN PMEUnit pvpme_unit; // polarization virial
TINKER_EXTERN PMEUnit dpme_unit;  // dispersion

TINKER_EXTERN real (*cmp)[10];
TINKER_EXTERN real (*fmp)[10];
TINKER_EXTERN real (*cphi)[10];
TINKER_EXTERN real (*fphi)[20];

TINKER_EXTERN real (*fuind)[3];
TINKER_EXTERN real (*fuinp)[3];
TINKER_EXTERN real (*fdip_phi1)[10];
TINKER_EXTERN real (*fdip_phi2)[10];
TINKER_EXTERN real (*cphidp)[10];
TINKER_EXTERN real (*fphidp)[20];

TINKER_EXTERN virial_buffer vir_m;
}
