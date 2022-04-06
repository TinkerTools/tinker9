#pragma once

namespace tinker {
/// \ingroup pme
/// \{
void pmeStreamStartRecord(bool usePmeStream);
void pmeStreamStartWait(bool usePmeStream);
void pmeStreamFinishRecord(bool usePmeStream);
void pmeStreamFinishWait(bool usePmeStream);
/// \}
}
