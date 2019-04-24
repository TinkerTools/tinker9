#ifndef TINKER_MOD_SOLUTE_HH_
#define TINKER_MOD_SOLUTE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace solute {
extern double& doffset;
extern double& p1;
extern double& p2;
extern double& p3;
extern double& p4;
extern double& p5;
extern double*& rsolv;
extern double*& asolv;
extern double*& rborn;
extern double*& drb;
extern double*& drbp;
extern double*& drobc;
extern double*& gpol;
extern double*& shct;
extern double*& aobc;
extern double*& bobc;
extern double*& gobc;
extern double*& vsolv;
extern double*& wace;
extern double*& s2ace;
extern double*& uace;
extern char (&solvtyp)[8];
extern char (&borntyp)[8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(solute, doffset);
extern "C" double m_tinker_mod(solute, p1);
extern "C" double m_tinker_mod(solute, p2);
extern "C" double m_tinker_mod(solute, p3);
extern "C" double m_tinker_mod(solute, p4);
extern "C" double m_tinker_mod(solute, p5);
extern "C" double* m_tinker_mod(solute, rsolv);
extern "C" double* m_tinker_mod(solute, asolv);
extern "C" double* m_tinker_mod(solute, rborn);
extern "C" double* m_tinker_mod(solute, drb);
extern "C" double* m_tinker_mod(solute, drbp);
extern "C" double* m_tinker_mod(solute, drobc);
extern "C" double* m_tinker_mod(solute, gpol);
extern "C" double* m_tinker_mod(solute, shct);
extern "C" double* m_tinker_mod(solute, aobc);
extern "C" double* m_tinker_mod(solute, bobc);
extern "C" double* m_tinker_mod(solute, gobc);
extern "C" double* m_tinker_mod(solute, vsolv);
extern "C" double* m_tinker_mod(solute, wace);
extern "C" double* m_tinker_mod(solute, s2ace);
extern "C" double* m_tinker_mod(solute, uace);
extern "C" char m_tinker_mod(solute, solvtyp)[8];
extern "C" char m_tinker_mod(solute, borntyp)[8];

double& doffset = m_tinker_mod(solute, doffset);
double& p1 = m_tinker_mod(solute, p1);
double& p2 = m_tinker_mod(solute, p2);
double& p3 = m_tinker_mod(solute, p3);
double& p4 = m_tinker_mod(solute, p4);
double& p5 = m_tinker_mod(solute, p5);
double*& rsolv = m_tinker_mod(solute, rsolv);
double*& asolv = m_tinker_mod(solute, asolv);
double*& rborn = m_tinker_mod(solute, rborn);
double*& drb = m_tinker_mod(solute, drb);
double*& drbp = m_tinker_mod(solute, drbp);
double*& drobc = m_tinker_mod(solute, drobc);
double*& gpol = m_tinker_mod(solute, gpol);
double*& shct = m_tinker_mod(solute, shct);
double*& aobc = m_tinker_mod(solute, aobc);
double*& bobc = m_tinker_mod(solute, bobc);
double*& gobc = m_tinker_mod(solute, gobc);
double*& vsolv = m_tinker_mod(solute, vsolv);
double*& wace = m_tinker_mod(solute, wace);
double*& s2ace = m_tinker_mod(solute, s2ace);
double*& uace = m_tinker_mod(solute, uace);
char (&solvtyp)[8] = m_tinker_mod(solute, solvtyp);
char (&borntyp)[8] = m_tinker_mod(solute, borntyp);
#endif

} TINKER_NAMESPACE_END

#endif
