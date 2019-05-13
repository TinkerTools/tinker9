#include "files.h"

#define FILE_DECL_(symb) const char* symb = current::symb

TINKER_NAMESPACE_BEGIN
namespace test {
namespace current = tinker_6fe8e913fe4da3d46849d10248ad2a4872b4da93;

// prm files
#include "file/water03.hh"
FILE_DECL_(water03_prm);
#include "file/amoeba09.hh"
FILE_DECL_(amoeba09_prm);

// xyz files
#include "file/watersmall.hh"
FILE_DECL_(watersmall_xyz);

// NaCl ion pair test set
#include "file/nacl.hh"
}
TINKER_NAMESPACE_END

#undef FILE_DECL_
