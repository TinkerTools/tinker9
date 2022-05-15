#pragma once
#include "ff/precision.h"

namespace tinker {
void dfieldAplus(real (*field)[3]);

void ufieldAplus(const real (*uind)[3], real (*field)[3]);

void diagPrecond3(const real (*rsd)[3], real (*zrsd)[3]);

void sparsePrecondBuild3();
void sparsePrecondApply3(const real (*rsd)[3], real (*zrsd)[3]);
}
