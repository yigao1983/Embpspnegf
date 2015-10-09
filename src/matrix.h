#ifndef matrix_h
#define matrix_h

#include "precision.h"

typedef struct mat1d_dreal
{
  int   ndim;
  dreal *addr;
  void (*alloc)(int, struct mat1d_dreal *);
  void (*deall)(struct mat1d_dreal *);
  void (*print)(const char *, struct mat1d_dreal *);
  void (*reset)(struct mat1d_dreal *);
} mat1d_dreal;

typedef struct mat1d_dcmplx
{
  int    ndim;
  dcmplx *addr;
  void (*alloc)(int, struct mat1d_dcmplx *);
  void (*deall)(struct mat1d_dcmplx *);
  void (*print)(const char *, struct mat1d_dcmplx *);
  void (*reset)(struct mat1d_dcmplx *);
} mat1d_dcmplx;

typedef struct mat2d_dreal
{
  int   nrow, ncol;
  dreal *addr;
  dreal **ptr;
  void (*alloc)(int, int, struct mat2d_dreal *);
  void (*deall)(struct mat2d_dreal *);
  void (*print)(const char *, struct mat2d_dreal *);
  void (*reset)(struct mat2d_dreal *);
} mat2d_dreal;

typedef struct mat2d_dcmplx
{
  int    nrow, ncol;
  dcmplx *addr;
  dcmplx **ptr;
  void (*alloc)(int, int, struct mat2d_dcmplx *);
  void (*deall)(struct mat2d_dcmplx *);
  void (*print)(const char *, struct mat2d_dcmplx *);
  void (*reset)(struct mat2d_dcmplx *);
} mat2d_dcmplx;

mat1d_dreal *mat1d_new_dreal(void);
mat1d_dcmplx *mat1d_new_dcmplx(void);

void mat1d_alloc_dreal(int, mat1d_dreal *);
void mat1d_alloc_dcmplx(int, mat1d_dcmplx *);

void mat1d_deall_dreal(mat1d_dreal *);
void mat1d_deall_dcmplx(mat1d_dcmplx *);

void mat1d_print_dreal(const char *, mat1d_dreal *);
void mat1d_print_dcmplx(const char *, mat1d_dcmplx *);

void mat1d_reset_dreal(mat1d_dreal *);
void mat1d_reset_dcmplx(mat1d_dcmplx *);

void mat1d_copy_dreal(mat1d_dreal *, mat1d_dreal *);
void mat1d_copy_dcmplx(mat1d_dcmplx *, mat1d_dcmplx *);

mat2d_dreal *mat2d_new_dreal(void);
mat2d_dcmplx *mat2d_new_dcmplx(void);

void mat2d_alloc_dreal(int, int, mat2d_dreal *);
void mat2d_alloc_dcmplx(int, int, mat2d_dcmplx *);

void mat2d_deall_dreal(mat2d_dreal *);
void mat2d_deall_dcmplx(mat2d_dcmplx *);

void mat2d_print_dreal(const char *, mat2d_dreal *);
void mat2d_print_dcmplx(const char *, mat2d_dcmplx *);

void mat2d_reset_dreal(mat2d_dreal *);
void mat2d_reset_dcmplx(mat2d_dcmplx *);

void mat2d_copy_dreal(mat2d_dreal *, mat2d_dreal *);
void mat2d_copy_dcmplx(mat2d_dcmplx *, mat2d_dcmplx *);

// For debug
//void mat2d_check(void);

#define mat1d_del_dreal(A) { mat1d_deall_dreal(A); freeup(A); }
#define mat1d_del_dcmplx(A) { mat1d_deall_dcmplx(A); freeup(A); }
#define mat2d_del_dreal(A) { mat2d_deall_dreal(A); freeup(A); }
#define mat2d_del_dcmplx(A) { mat2d_deall_dcmplx(A); freeup(A); }

#endif
