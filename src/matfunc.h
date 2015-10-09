#ifndef matfunc_h
#define matfunc_h

#include "matrix.h"

typedef struct matfunc_dreal
{
  int         ngrd, ntag, nrow, ncol;
  int         *tag;
  dreal       *xval;
  mat2d_dreal **ptr;
  void (*alloc)(int, int, int, struct matfunc_dreal *);
  void (*deall)(struct matfunc_dreal *);
} matfunc_dreal;

typedef struct matfunc_dcmplx
{
  int          ngrd, ntag, nrow, ncol;
  int          *tag;
  dreal        *xval;
  mat2d_dcmplx **ptr;
  void (*alloc)(int, int, int, struct matfunc_dcmplx *);
  void (*deall)(struct matfunc_dcmplx *);
} matfunc_dcmplx;

matfunc_dreal *matfunc_new_dreal(void);
void matfunc_alloc_dreal(int, int, int, matfunc_dreal *);
void matfunc_deall_dreal(matfunc_dreal *);
void matfunc_copy_dreal(matfunc_dreal *, matfunc_dreal *);
void matfunc_bisec_dreal(dreal, matfunc_dreal *, mat2d_dreal *);
void matfunc_intrp_tag_dreal(matfunc_dreal *);

matfunc_dcmplx *matfunc_new_dcmplx(void);
void matfunc_alloc_dcmplx(int, int, int, matfunc_dcmplx *);
void matfunc_deall_dcmplx(matfunc_dcmplx *);
void matfunc_copy_dcmplx(matfunc_dcmplx *, matfunc_dcmplx *);
void matfunc_bisec_dcmplx(dreal, matfunc_dcmplx *, mat2d_dcmplx *);
void matfunc_intrp_tag_dcmplx(matfunc_dcmplx *);

// For debug
//void matfunc_check_dreal(void);
//void matfunc_check_dcmplx(void);

#define matfunc_del_dreal(A) { matfunc_deall_dreal(A); freeup(A); }
#define matfunc_del_dcmplx(A) { matfunc_deall_dcmplx(A); freeup(A); }

#endif
