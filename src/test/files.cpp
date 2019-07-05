#include "files.h"

#define FILE_DEF_(symb) const char* symb = current::symb

TINKER_NAMESPACE_BEGIN
namespace test {
namespace current = tinker_6fe8e913fe4da3d46849d10248ad2a4872b4da93;

// prm files
#include "file/water03.hh"
FILE_DEF_(water03_prm);
#include "file/amoeba09.hh"
FILE_DEF_(amoeba09_prm);
#include "file/amoebabio09.hh"
FILE_DEF_(amoebabio09_prm);
#include "file/amoebapro13.hh"
FILE_DEF_(amoebapro13_prm);

// xyz files
#include "file/watersmall.hh"
FILE_DEF_(watersmall_xyz);

// NaCl ion pair test set
#include "file/nacl.hh"

// CLN025 test set
#include "file/cln025.hh"

// local_frame test set
#include "file/local_frame.hh"

// crambin3 test set
#include "file/crambin3.hh"

// trpcage test set
#include "file/trpcage.hh"

// ten water molecules test set
#include "file/h2o10.hh"

// six argon atoms test set
#include "file/ar6.hh"
}
TINKER_NAMESPACE_END

#undef FILE_DEF_
