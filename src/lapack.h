#ifndef lapack_h
#define lapack_h

#include "precision.h"

extern void dgemm_(char *, char *, int *, int *, int *,
		   dreal *, dreal *, int *, dreal *, int *,
		   dreal *, dreal *, int *);
extern void dgetrf_(int *, int *, dreal *, int *, int *, int *);
extern void dgetri_(int *, dreal *, int *, int *,
                    dreal *, int *, int *);
extern void zgemm_(char *, char *, int *, int *, int *,
		   dcmplx *, dcmplx *, int *,
		   dcmplx *, int *, dcmplx *,
		   dcmplx *, int *);
extern void zgetrf_(int *, int *, dcmplx *, int *, int *, int *);
extern void zgetri_(int *, dcmplx *, int *, int *,
		    dcmplx *, int *, int *);

void lapack_dgemm(int, int, int, char, char,
		  dreal, dreal *, dreal *, dreal, dreal *);
void lapack_dgetrf(int, int, dreal *, int *);
void lapack_dgetri(int, dreal *, int *);
void lapack_dinvs(int, dreal *);

void lapack_zgemm(int, int, int, char, char,
		  dcmplx, dcmplx *, dcmplx *, dcmplx, dcmplx *);
void lapack_zgetrf(int, int, dcmplx *, int *);
void lapack_zgetri(int, dcmplx *, int *);
void lapack_zinvs(int, dcmplx *);

void lapack_check(void);

#endif
