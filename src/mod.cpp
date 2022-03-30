#define TINKER_EXTERN_DEFINITION_FILE 1

#include "ff/atom.h"
#include "ff/box.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/molecule.h"
#include "ff/nblist.h"
#include "ff/pme.h"
#include "ff/rattle.h"
#include "ff/spatial.h"

#include "ff/amoeba/elecamoeba.h"

#include "ff/hippo/edisp.h"
#include "ff/hippo/elechippo.h"
#include "ff/hippo/erepel.h"

#include "ff/pchg/echarge.h"
#include "ff/pchg/echglj.h"
#include "ff/pchg/evalence.h"
#include "ff/pchg/evdw.h"

#include "md/intg.h"
#include "md/pq.h"
#include "md/pt.h"

#include "tool/cudalib.h"
#include "tool/energybuffer.h"
#include "tool/gpucard.h"
#include "tool/rcman.h"

#include "accasync.h"
#include "platform.h"
