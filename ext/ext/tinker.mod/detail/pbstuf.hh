#ifndef TINKER_MOD_PBSTUF_HH_
#define TINKER_MOD_PBSTUF_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace pbstuf {
const int maxion = 10;
extern int& ionn;
extern int (&dime)[3];
extern int (&ionq)[maxion];
extern double& pbe;
extern double& pdie;
extern double& sdie;
extern double& srad;
extern double& swin;
extern double& sdens;
extern double& smin;
extern double (&grid)[3];
extern double (&gcent)[3];
extern double (&cgrid)[3];
extern double (&cgcent)[3];
extern double (&fgrid)[3];
extern double (&fgcent)[3];
extern double (&ionr)[maxion];
extern double (&ionc)[maxion];
extern double*& apbe;
extern double*& pbr;
extern double*& pbep;
extern double*& pbfp;
extern double*& pbtp;
extern double*& pbeuind;
extern double*& pbeuinp;
extern char (&pbtyp)[20];
extern char (&pbsoln)[20];
extern char (&bcfl)[20];
extern char (&chgm)[20];
extern char (&srfm)[20];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(pbstuf, ionn);
extern "C" int m_tinker_mod(pbstuf, dime)[3];
extern "C" int m_tinker_mod(pbstuf, ionq)[maxion];
extern "C" double m_tinker_mod(pbstuf, pbe);
extern "C" double m_tinker_mod(pbstuf, pdie);
extern "C" double m_tinker_mod(pbstuf, sdie);
extern "C" double m_tinker_mod(pbstuf, srad);
extern "C" double m_tinker_mod(pbstuf, swin);
extern "C" double m_tinker_mod(pbstuf, sdens);
extern "C" double m_tinker_mod(pbstuf, smin);
extern "C" double m_tinker_mod(pbstuf, grid)[3];
extern "C" double m_tinker_mod(pbstuf, gcent)[3];
extern "C" double m_tinker_mod(pbstuf, cgrid)[3];
extern "C" double m_tinker_mod(pbstuf, cgcent)[3];
extern "C" double m_tinker_mod(pbstuf, fgrid)[3];
extern "C" double m_tinker_mod(pbstuf, fgcent)[3];
extern "C" double m_tinker_mod(pbstuf, ionr)[maxion];
extern "C" double m_tinker_mod(pbstuf, ionc)[maxion];
extern "C" double* m_tinker_mod(pbstuf, apbe);
extern "C" double* m_tinker_mod(pbstuf, pbr);
extern "C" double* m_tinker_mod(pbstuf, pbep);
extern "C" double* m_tinker_mod(pbstuf, pbfp);
extern "C" double* m_tinker_mod(pbstuf, pbtp);
extern "C" double* m_tinker_mod(pbstuf, pbeuind);
extern "C" double* m_tinker_mod(pbstuf, pbeuinp);
extern "C" char m_tinker_mod(pbstuf, pbtyp)[20];
extern "C" char m_tinker_mod(pbstuf, pbsoln)[20];
extern "C" char m_tinker_mod(pbstuf, bcfl)[20];
extern "C" char m_tinker_mod(pbstuf, chgm)[20];
extern "C" char m_tinker_mod(pbstuf, srfm)[20];

int& ionn = m_tinker_mod(pbstuf, ionn);
int (&dime)[3] = m_tinker_mod(pbstuf, dime);
int (&ionq)[maxion] = m_tinker_mod(pbstuf, ionq);
double& pbe = m_tinker_mod(pbstuf, pbe);
double& pdie = m_tinker_mod(pbstuf, pdie);
double& sdie = m_tinker_mod(pbstuf, sdie);
double& srad = m_tinker_mod(pbstuf, srad);
double& swin = m_tinker_mod(pbstuf, swin);
double& sdens = m_tinker_mod(pbstuf, sdens);
double& smin = m_tinker_mod(pbstuf, smin);
double (&grid)[3] = m_tinker_mod(pbstuf, grid);
double (&gcent)[3] = m_tinker_mod(pbstuf, gcent);
double (&cgrid)[3] = m_tinker_mod(pbstuf, cgrid);
double (&cgcent)[3] = m_tinker_mod(pbstuf, cgcent);
double (&fgrid)[3] = m_tinker_mod(pbstuf, fgrid);
double (&fgcent)[3] = m_tinker_mod(pbstuf, fgcent);
double (&ionr)[maxion] = m_tinker_mod(pbstuf, ionr);
double (&ionc)[maxion] = m_tinker_mod(pbstuf, ionc);
double*& apbe = m_tinker_mod(pbstuf, apbe);
double*& pbr = m_tinker_mod(pbstuf, pbr);
double*& pbep = m_tinker_mod(pbstuf, pbep);
double*& pbfp = m_tinker_mod(pbstuf, pbfp);
double*& pbtp = m_tinker_mod(pbstuf, pbtp);
double*& pbeuind = m_tinker_mod(pbstuf, pbeuind);
double*& pbeuinp = m_tinker_mod(pbstuf, pbeuinp);
char (&pbtyp)[20] = m_tinker_mod(pbstuf, pbtyp);
char (&pbsoln)[20] = m_tinker_mod(pbstuf, pbsoln);
char (&bcfl)[20] = m_tinker_mod(pbstuf, bcfl);
char (&chgm)[20] = m_tinker_mod(pbstuf, chgm);
char (&srfm)[20] = m_tinker_mod(pbstuf, srfm);
#endif

} TINKER_NAMESPACE_END

#endif
