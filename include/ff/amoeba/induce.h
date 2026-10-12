#pragma once
#include "ff/modamoeba.h"
#include "ff/precision.h"
#include "tool/rcman.h"

#include <vector>

namespace tinker {
/// \ingroup polar
/// \{
// electrostatic field due to permanent multipoles
void dfield(real (*field)[3], real (*fieldp)[3]);
void dfieldNonEwald(real (*field)[3], real (*fieldp)[3]);
void dfieldEwald(real (*field)[3], real (*fieldp)[3]);
void dfieldEwaldRecipSelfP1(real (*field)[3]);

// mutual electrostatic field due to induced dipole moments
// -Tu operator
void ufield(const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufieldNonEwald(const real (*uind)[3],
                    const real (*uinp)[3], //
                    real (*field)[3],
                    real (*fieldp)[3]);
void ufieldEwald(const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);

void diagPrecond(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3]);

void sparsePrecondBuild();
void sparsePrecondApply(const real (*rsd)[3],
                        const real (*rsdp)[3], //
                        real (*zrsd)[3],
                        real (*zrsdp)[3]);

void ulspredSave(const real (*uind)[3], const real (*uinp)[3]);

/// The induced dipole predictor a run uses: the choice the Fortran reader made
/// (predict.f), except for the dynamics of free energy perturbation and of the
/// adaptive lambda bias methods, where the reader turns prediction off and the
/// POLAR-PREDICT keyword is read again here.
UPred polarPredictor();
void ulspredSum(real (*uind)[3], real (*uinp)[3]);

/// A copy of the induced dipoles, the predictor history and, with exchange
/// polarization, its matrices, so that an energy evaluated at a trial geometry
/// which is then rejected leaves no trace.
class InducedSnapshot
{
public:
   InducedSnapshot();
   ~InducedSnapshot();
   InducedSnapshot(const InducedSnapshot&) = delete;
   InducedSnapshot& operator=(const InducedSnapshot&) = delete;
   /// \c history is false when the solves to come cannot write the predictor
   /// history, which leaves its slot out of the copy.
   void save(bool history = true);
   void restore();

private:
   template <class T>
   struct Copy
   {
      T** src;
      T* copy;
   };
   std::vector<Copy<real[3]>> m_vec;
   std::vector<Copy<real[3][3]>> m_mat;
   Copy<real[3]> m_ud, m_up; // the one predictor slot a trial overwrites
   bool m_slot;              // whether the last save() copied it
   int m_nualt;
};

/// An energy evaluated at a state the dynamics is not in. When it ends the
/// induced dipoles, and whatever its solves changed of the predictor history,
/// are put back, unless accept() was called. Trials do not nest.
class InducedTrial
{
public:
   /// What the solves of the trial may do with the predictor history.
   enum class History
   {
      USE,   ///< Read and write it as the dynamics does: a trial geometry close
             ///< to the current one, such as a Monte Carlo volume move.
      BYPASS ///< Neither read nor write it: they start from the direct guess,
             ///< their result is not saved, and one that fails to converge is
             ///< not retried by wiggling the atoms. For the trial lambdas of
             ///< free energy perturbation and of a Monte Carlo lambda move.
   };

   explicit InducedTrial(History);
   ~InducedTrial();
   InducedTrial(const InducedTrial&) = delete;
   InducedTrial& operator=(const InducedTrial&) = delete;
   /// Keeps what the trial left behind, as for an accepted move.
   void accept();

private:
   bool m_saved, m_accepted;
};

/// How an induced dipole solve starts.
enum class UGuess
{
   PREDICT, ///< from the predictor history (ulspredSum)
   DIRECT,  ///< from the direct dipoles
   ZERO     ///< from zero
};

/// The starting guess of the solve at hand: the predictor when it is on, its
/// history is full, and neither a wiggle retry nor a trial bypasses it; the
/// direct guess when the predictor had to be passed over or PCG-GUESS asks for
/// it; zero otherwise.
UGuess induceGuess();

void inducePrint(const real (*ud)[3]);
void induce(real (*uind)[3], real (*uinp)[3]);

// recovery from unconverged induced dipoles in MD (keyword INDUCE-WIGGLE)
void induceWiggleData(RcOp);
/// For a solver that ran out of iterations. Under INDUCE-WIGGLE in MD the
/// failure is reported and recorded for induceWiggleRetry(), and the solver
/// returns; otherwise, or within a trial that bypasses the history, this
/// throws. \c rsdp may be null.
void induceNotConverged(const real (*rsd)[3], const real (*rsdp)[3]);
bool induceFailed();          // a solver failed in the current energy evaluation
bool induceHistoryBypassed(); // within an InducedTrial with History::BYPASS
bool induceWiggleRetry();     // wiggles atoms and returns true if energy must be recomputed
/// \}
}

//          | h        | cpp        | acc            | cu
// field    | induce.h | field.cpp  | acc/field*.cpp | field.cu
// pcg      | --       | --         | acc/induce.cpp | pcg.cu
// precond  | induce.h | induce.cpp | acc/induce.cpp | precond.cu
// upredict | induce.h | induce.cpp | acc/induce.cpp | --
