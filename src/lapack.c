#include <stdio.h>
#include <stdlib.h>
#include "dbg.h"
#include "sys.h"
#include "matrix.h"
#include "lapack.h"

// Interface to lapack routine dgemm
// mm: nrow of A
// nn: ncol of B
// kk: ncol and nrow of C
void lapack_dgemm(int mm, int nn, int kk,
                  char transa, char transb,
                  dreal alpha, dreal *AA, dreal *BB,
                  dreal beta, dreal *CC)
{
  int lda, ldb, ldc;
  
  if(transa == 'N' || transa == 'n') {
    lda = 1 > mm ? 1 : mm;
  } else {
    lda = 1 > kk ? 1 : kk;
  }
  
  if(transb == 'N' || transb == 'n') {
    ldb = 1 > kk ? 1 : kk;
  } else {
    ldb = 1 > nn ? 1 : nn;
  }
  
  ldc = 1 > mm ? 1 : mm;
  
  dgemm_(&transa, &transb, &mm, &nn, &kk, &alpha, AA, &lda, BB, &ldb, &beta, CC, &ldc); 
}

// Interface to lapack routine dgetrf
// mm: nrow of A
// nn: ncol of A
void lapack_dgetrf(int mm, int nn, dreal *AA, int *ipiv)
{
  int lda, info;
  
  lda = 1 > mm ? 1 : mm;
  
  dgetrf_(&mm, &nn, AA, &lda, ipiv, &info);
  check(info == 0, "Failed dgetrf, info = %d", info);
  
  return;
  
 error:
  abort();
}

// Interface to lapack routine dgetri
// nn: nrow and ncol of A
void lapack_dgetri(int nn, dreal *AA, int *ipiv)
{
  int   lda, lwork, info;
  dreal *work = NULL;
  
  lda = 1 > nn ? 1 : nn;
  lwork = lda;
  
  work = (dreal *) calloc(lwork, sizeof(dreal));
  check_mem(work, "work");
  
  dgetri_(&nn, AA, &lda, ipiv, work, &lwork, &info);
  check(info == 0, "Failed dgetri, info = %d", info);
  
  freeup(work);
  
  return;
  
 error:
  if(work) freeup(work);
  abort();
}

// Inverse of matrix AA
void lapack_dinvs(int nn, dreal *AA)
{
  int *ipiv = NULL;
  
  ipiv = (int *) calloc(nn, sizeof(int));
  check_mem(ipiv, "ipiv");
  
  lapack_dgetrf(nn, nn, AA, ipiv);
  lapack_dgetri(nn, AA, ipiv);
  
  freeup(ipiv);
  
  return;
  
 error:
  if(ipiv) freeup(ipiv);
  abort();
}

// Interface to lapack routine zgemm
// mm: nrow of A
// nn: ncol of B
// kk: ncol and nrow of C
void lapack_zgemm(int mm, int nn, int kk,
                  char transa, char transb,
                  dcmplx alpha, dcmplx *AA, dcmplx *BB,
                  dcmplx beta, dcmplx *CC)
{
  int lda, ldb, ldc;
  
  if(transa == 'N' || transa == 'n') {
    lda = 1 > mm ? 1 : mm;
  } else {
    lda = 1 > kk ? 1 : kk;
  }
  
  if(transb == 'N' || transb == 'n') {
    ldb = 1 > kk ? 1 : kk;
  } else {
    ldb = 1 > nn ? 1 : nn;
  }
  
  ldc = 1 > mm ? 1 : mm;
  
  zgemm_(&transa, &transb, &mm, &nn, &kk, &alpha, AA, &lda, BB, &ldb, &beta, CC, &ldc); 
}

// Interface to lapack routine zgetrf
// mm: nrow of A
// nn: ncol of A
void lapack_zgetrf(int mm, int nn, dcmplx *AA, int *ipiv)
{
  int lda, info;
  
  lda = 1 > mm ? 1 : mm;
  
  zgetrf_(&mm, &nn, AA, &lda, ipiv, &info);
  check(info == 0, "Failed zgetrf, info = %d", info);
  
  return;
  
 error:
  abort();
}

// Interface to lapack routine zgetri
// nn: nrow and ncol of A
void lapack_zgetri(int nn, dcmplx *AA, int *ipiv)
{
  int     lda, lwork, info;
  dcmplx  *work = NULL;
  
  lda = 1 > nn ? 1 : nn;
  lwork = lda;
  
  work = (dcmplx *) calloc(lwork, sizeof(dcmplx));
  check_mem(work, "work");
  
  zgetri_(&nn, AA, &lda, ipiv, work, &lwork, &info);
  check(info == 0, "Failed zgetri, info = %d", info);
  
  freeup(work);
  
  return;
  
 error:
  if(work) freeup(work);
  abort();
}

// Inverse of matrix AA
void lapack_zinvs(int nn, dcmplx *AA)
{
  int *ipiv = NULL;
  
  ipiv = (int *) calloc(nn, sizeof(int));
  check_mem(ipiv, "ipiv");
  
  lapack_zgetrf(nn, nn, AA, ipiv);
  lapack_zgetri(nn, AA, ipiv);
  
  freeup(ipiv);
  
  return;
  
 error:
  if(ipiv) freeup(ipiv);
  abort();
}
