#ifndef TINKER_MOD_PDB_HH_
#define TINKER_MOD_PDB_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace pdb {
extern int& npdb;
extern int& nres;
extern int*& resnum;
extern int*& resatm;
extern int*& npdb12;
extern int*& ipdb12;
extern int*& pdblist;
extern double*& xpdb;
extern double*& ypdb;
extern double*& zpdb;
extern char (&altsym)[1];
extern char (*&pdbres)[3];
extern char (*&pdbatm)[4];
extern char (*&pdbtyp)[6];
extern char (&chnsym)[20];
extern char (&instyp)[20];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(pdb, npdb);
extern "C" int m_tinker_mod(pdb, nres);
extern "C" int* m_tinker_mod(pdb, resnum);
extern "C" int* m_tinker_mod(pdb, resatm);
extern "C" int* m_tinker_mod(pdb, npdb12);
extern "C" int* m_tinker_mod(pdb, ipdb12);
extern "C" int* m_tinker_mod(pdb, pdblist);
extern "C" double* m_tinker_mod(pdb, xpdb);
extern "C" double* m_tinker_mod(pdb, ypdb);
extern "C" double* m_tinker_mod(pdb, zpdb);
extern "C" char m_tinker_mod(pdb, altsym)[1];
extern "C" char (*m_tinker_mod(pdb, pdbres))[3];
extern "C" char (*m_tinker_mod(pdb, pdbatm))[4];
extern "C" char (*m_tinker_mod(pdb, pdbtyp))[6];
extern "C" char m_tinker_mod(pdb, chnsym)[20];
extern "C" char m_tinker_mod(pdb, instyp)[20];

int& npdb = m_tinker_mod(pdb, npdb);
int& nres = m_tinker_mod(pdb, nres);
int*& resnum = m_tinker_mod(pdb, resnum);
int*& resatm = m_tinker_mod(pdb, resatm);
int*& npdb12 = m_tinker_mod(pdb, npdb12);
int*& ipdb12 = m_tinker_mod(pdb, ipdb12);
int*& pdblist = m_tinker_mod(pdb, pdblist);
double*& xpdb = m_tinker_mod(pdb, xpdb);
double*& ypdb = m_tinker_mod(pdb, ypdb);
double*& zpdb = m_tinker_mod(pdb, zpdb);
char (&altsym)[1] = m_tinker_mod(pdb, altsym);
char (*&pdbres)[3] = m_tinker_mod(pdb, pdbres);
char (*&pdbatm)[4] = m_tinker_mod(pdb, pdbatm);
char (*&pdbtyp)[6] = m_tinker_mod(pdb, pdbtyp);
char (&chnsym)[20] = m_tinker_mod(pdb, chnsym);
char (&instyp)[20] = m_tinker_mod(pdb, instyp);
#endif

} TINKER_NAMESPACE_END

#endif
