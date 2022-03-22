#pragma once
#include "tool/rcman.h"

namespace tinker {
void evalence(int vers);

//

enum class ebond_t
{
   harmonic,
   morse
};
void ebondData(RcOp);
void ebond(int vers);
void ebond_acc(int);

//

enum class eangle_t : int
{
   in_plane,
   harmonic,
   linear,
   fourier
};
enum class eopbend_t
{
   w_d_c,
   allinger
};
void eangleData(RcOp);
void eangle(int vers);
void eangle_acc(int);

//

void estrbndData(RcOp);
void estrbnd(int vers);
void estrbnd_acc(int);

//

void eureyData(RcOp);
void eurey(int vers);
void eurey_acc(int);

//

void eopbendData(RcOp);
void eopbend(int vers);
void eopbend_acc(int);

//

void eimpropData(RcOp);
void eimprop(int vers);
void eimprop_acc(int);

//

void eimptorData(RcOp);
void eimptor(int vers);
void eimptor_acc(int);

//

void etorsData(RcOp);
void etors(int vers);
void etors_acc(int);

//

void epitorsData(RcOp);
void epitors(int vers);
void epitors_acc(int);

//

void estrtorData(RcOp);
void estrtor(int vers);
void estrtor_acc(int);

//

void eangtorData(RcOp);
void eangtor(int vers);
void eangtor_acc(int);

//

void etortorData(RcOp);
void etortor(int vers);
void etortor_acc(int);

//

void egeomData(RcOp);
void egeom(int vers);
void egeom_acc(int);
}

#include "mod/evalence.h"
