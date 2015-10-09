#ifndef numeric_h
#define numeric_h

#include "precision.h"

dreal numeric_maxval_dreal(int, dreal *);
dreal numeric_minval_dreal(int, dreal *);
dreal numeric_quadrat_dreal(int, dreal *, dreal *);
dcmplx numeric_quadrat_dcmplx(int, dreal *, dcmplx *);
int numeric_ipow(int, int);
int numeric_ispow2(int);
void numeric_kronig(int, dcmplx *, dcmplx *);
void numeric_convolut(int, dreal, int, dcmplx *, int, dcmplx *, dcmplx *);
void numeric_mirror(int, dcmplx *, dcmplx *);
void numeric_freq2time(int, dreal, dcmplx *, dcmplx *);
void numeric_time2freq(int, dreal, dcmplx *, dcmplx *);

#endif
