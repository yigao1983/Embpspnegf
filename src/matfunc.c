#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dbg.h"
#include "sys.h"
#include "constants.h"
#include "matfunc.h"

matfunc_dreal *matfunc_new_dreal(void)
{
  matfunc_dreal *matfunc = NULL;
  
  matfunc = (matfunc_dreal *) calloc(1, sizeof(matfunc_dreal));
  check_mem(matfunc, "matfunc");
  
  matfunc->alloc = matfunc_alloc_dreal;
  matfunc->deall = matfunc_deall_dreal;
  
  return matfunc;
  
 error:
  if(matfunc) freeup(matfunc);
  return NULL;
}

void matfunc_alloc_dreal(int ngrd, int nrow, int ncol, matfunc_dreal *matfunc)
{
  int         igrd;
  int         *tag  = NULL;
  dreal       *xval = NULL;
  mat2d_dreal **ptr = NULL;
  
  check(ngrd >= 0, "Invalid ngrd: %5d", ngrd);
  check(nrow >= 0, "Invalid nrow: %5d", nrow);
  check(ncol >= 0, "Invalid ncol: %5d", ncol);
  
  tag = (int *) calloc(ngrd, sizeof(int));
  check_mem(tag, "tag");
  xval = (dreal *) calloc(ngrd, sizeof(dreal));
  check_mem(xval, "xval");
  ptr = (mat2d_dreal **) calloc(ngrd, sizeof(mat2d_dreal *));
  check_mem(ptr, "ptr");
  
  for(igrd = 0; igrd < ngrd; ++igrd) {
    ptr[igrd] = mat2d_new_dreal();
    ptr[igrd]->alloc(nrow, ncol, ptr[igrd]);
  }
  
  matfunc->ngrd = ngrd;
  matfunc->nrow = nrow;
  matfunc->ncol = ncol;
  matfunc->tag  = tag;
  matfunc->xval = xval;
  matfunc->ptr  = ptr;
  
  return;
  
 error:
  if(tag)  freeup(tag);
  if(xval) freeup(xval);
  if(ptr)  freeup(ptr);
  abort();
}

void matfunc_deall_dreal(matfunc_dreal *matfunc)
{
  int ngrd, igrd;
  
  if(matfunc) {
    ngrd = matfunc->ngrd;
    
    for(igrd = ngrd - 1; igrd > -1; --igrd) {
      mat2d_del_dreal(matfunc->ptr[igrd]);
    }
    
    if(matfunc->ptr)  freeup(matfunc->ptr);
    if(matfunc->xval) freeup(matfunc->xval);
    if(matfunc->tag)  freeup(matfunc->tag);
  }
}

void matfunc_copy_dreal(matfunc_dreal *matfunc_src, matfunc_dreal *matfunc_des)
{
  int    ngrd_src, ngrd_des, igrd, nrow_src, nrow_des, ncol_src, ncol_des;
  size_t nlen;
  
  check_mem(matfunc_src, "matfunc_src");
  check_mem(matfunc_des, "matfunc_des");
  ngrd_src = matfunc_src->ngrd; ngrd_des = matfunc_des->ngrd;
  nrow_src = matfunc_src->nrow; nrow_des = matfunc_des->nrow;
  ncol_src = matfunc_src->ncol; ncol_des = matfunc_des->ncol;
  check(ngrd_src == ngrd_des, "Inequivalent ngrd");
  check(nrow_src == nrow_des, "Inequivalent nrow");
  check(ncol_src == ncol_des, "Inequivalent ncol");
  
  matfunc_des->ntag = matfunc_src->ntag;
  
  nlen = ngrd_src * sizeof(int);
  memcpy(matfunc_des->tag, matfunc_src->tag, nlen);
  
  nlen = ngrd_src * sizeof(dreal);
  memcpy(matfunc_des->xval, matfunc_src->xval, nlen);
  
  for(igrd = 0; igrd < ngrd_src; ++igrd) {
    mat2d_copy_dreal(matfunc_src->ptr[igrd], matfunc_des->ptr[igrd]);
  }
  
  return;
  
 error:
  abort();
}

void matfunc_bisec_dreal(dreal xx, matfunc_dreal *matfunc, mat2d_dreal *mat)
{
  int         ngrd, igl, igr, igm, nrow, irow, ncol, icol;
  dreal       xl, xr, xm;
  mat2d_dreal *matl = NULL, *matr = NULL;
  
  ngrd = matfunc->ngrd;
  nrow = matfunc->nrow;
  ncol = matfunc->nrow;
  
  igl = 0;        xl = matfunc->xval[igl];
  igr = ngrd - 1; xr = matfunc->xval[igr];
  
  check((xx-xl) * (xx-xr) <= 0.0, "Invalid xx: %15.5e\n", xx);
  
  while(igr - igl > 1) {
    igm = (igl + igr) / 2; xm = matfunc->xval[igm];
    if((xx-xl) * (xx-xm) <= 0.0) {
      igr = igm; xr = xm;
    } else if((xx-xm) * (xx-xr) <= 0.0) {
      igl = igm; xl = xm;
    }
  }
  
  matl = matfunc->ptr[igl];
  matr = matfunc->ptr[igr];
  for(icol = 0; icol < ncol; ++icol) {
    for(irow = 0; irow < nrow; ++irow) {
      mat->ptr[icol][irow] = matl->ptr[icol][irow]
        + (xx - xl) / (xr - xl) *(matr->ptr[icol][irow] - matl->ptr[icol][irow]);
    }
  }
  
  return;
  
 error:
  abort();
}

void matfunc_intrp_tag_dreal(matfunc_dreal *matfunc)
{
  int         ngrd, ig, nrow, irow, ncol, icol, ntag, nlen;
  dreal       fv, fmax, df;
  dreal       fl, fr;
  mat2d_dreal *mat = NULL, *matl = NULL, *matr = NULL;
  
  ngrd = matfunc->ngrd;
  nrow = matfunc->nrow;
  ncol = matfunc->ncol;
  
  check(ngrd > 1, "Invalid ngrd: %5d", ngrd);
  
  nlen = ngrd * sizeof(int);
  memset(matfunc->tag, 0, nlen);
  matfunc->ntag = 0;
  
  fmax = 0.0;
  for(ig = 0; ig < ngrd; ++ig) {
    mat = matfunc->ptr[ig];
    for(icol = 0; icol < ncol; ++icol) {
      for(irow = 0; irow < nrow; ++irow) {
        fv = fabs(mat->ptr[icol][irow]);
        if(fv > fmax) fmax = fv;
      }
    }
  }
  
  for(ig = 0, ntag = 0; ig < ngrd - 1; ++ig) {
    matl = matfunc->ptr[ig];
    matr = matfunc->ptr[ig+1];
    for(icol = 0; icol < ncol; ++icol) {
      for(irow = 0; irow < nrow; ++irow) {
        fl = matl->ptr[icol][irow];
        fr = matr->ptr[icol][irow];
        df = fr - fl;
        if(fabs(df) > RATIO * fmax) {
          matfunc->tag[ig] = 1;
        }
      }
    }
    if(matfunc->tag[ig] == 1) ++ntag;
  }
  matfunc->ntag = ntag;
  
  return;
  
 error:
  abort();
}

matfunc_dcmplx *matfunc_new_dcmplx(void)
{
  matfunc_dcmplx *matfunc = NULL;
  
  matfunc = (matfunc_dcmplx *) calloc(1, sizeof(matfunc_dcmplx));
  check_mem(matfunc, "matfunc");
  
  matfunc->alloc = matfunc_alloc_dcmplx;
  matfunc->deall = matfunc_deall_dcmplx;
  
  return matfunc;
  
 error:
  if(matfunc) freeup(matfunc);
  return NULL;
}

void matfunc_alloc_dcmplx(int ngrd, int nrow, int ncol, matfunc_dcmplx *matfunc)
{
  int          igrd;
  int          *tag  = NULL;
  dreal        *xval = NULL;
  mat2d_dcmplx **ptr = NULL;
  
  check(ngrd >= 0, "Invalid ngrd: %5d", ngrd);
  check(nrow >= 0, "Invalid nrow: %5d", nrow);
  check(ncol >= 0, "Invalid ncol: %5d", ncol);
  
  tag  = (int *) calloc(ngrd, sizeof(int));
  check_mem(tag, "tag");
  xval = (dreal *) calloc(ngrd, sizeof(dreal));
  check_mem(xval, "xval");
  ptr  = (mat2d_dcmplx **) calloc(ngrd, sizeof(mat2d_dcmplx *));
  check_mem(ptr, "ptr");
  
  for(igrd = 0; igrd < ngrd; ++igrd) {
    ptr[igrd] = mat2d_new_dcmplx();
    ptr[igrd]->alloc(nrow, ncol, ptr[igrd]);
  }
  
  matfunc->ngrd = ngrd;
  matfunc->nrow = nrow;
  matfunc->ncol = ncol;
  matfunc->tag  = tag;
  matfunc->xval = xval;
  matfunc->ptr  = ptr;
  
  return;
  
 error:
  if(tag)  freeup(tag);
  if(xval) freeup(xval);
  if(ptr)  freeup(ptr);
  abort();
}

void matfunc_deall_dcmplx(matfunc_dcmplx *matfunc)
{
  int ngrd, igrd;
  
  if(matfunc) {
    ngrd = matfunc->ngrd;
    
    for(igrd = ngrd - 1; igrd > -1; --igrd) {
      mat2d_del_dcmplx(matfunc->ptr[igrd]);
    }
    
    if(matfunc->ptr)  freeup(matfunc->ptr);
    if(matfunc->xval) freeup(matfunc->xval);
    if(matfunc->tag)  freeup(matfunc->tag);
  }
}

void matfunc_copy_dcmplx(matfunc_dcmplx *matfunc_src, matfunc_dcmplx *matfunc_des)
{
  int    ngrd_src, ngrd_des, igrd, nrow_src, nrow_des, ncol_src, ncol_des;
  size_t nlen;
  
  check_mem(matfunc_src, "matfunc_src");
  check_mem(matfunc_des, "matfunc_des");
  ngrd_src = matfunc_src->ngrd; ngrd_des = matfunc_des->ngrd;
  nrow_src = matfunc_src->nrow; nrow_des = matfunc_des->nrow;
  ncol_src = matfunc_src->ncol; ncol_des = matfunc_des->ncol;
  check(ngrd_src == ngrd_des, "Inequivalent ngrd");
  check(nrow_src == nrow_des, "Inequivalent nrow");
  check(ncol_src == ncol_des, "Inequivalent ncol");
  
  matfunc_des->ntag = matfunc_src->ntag;
  
  nlen = ngrd_src * sizeof(int);
  memcpy(matfunc_des->tag, matfunc_src->tag, nlen);
  
  nlen = ngrd_src * sizeof(dreal);
  memcpy(matfunc_des->xval, matfunc_src->xval, nlen);
  
  for(igrd = 0; igrd < ngrd_src; ++igrd) {
    mat2d_copy_dcmplx(matfunc_src->ptr[igrd], matfunc_des->ptr[igrd]);
  }
  
  return;
  
 error:
  abort();
}

void matfunc_bisec_dcmplx(dreal xx, matfunc_dcmplx *matfunc, mat2d_dcmplx *mat)
{
  int          ngrd, igl, igr, igm, nrow, irow, ncol, icol;
  dreal        xl, xr, xm;
  mat2d_dcmplx *matl = NULL, *matr = NULL;
  
  ngrd = matfunc->ngrd;
  nrow = matfunc->nrow;
  ncol = matfunc->nrow;
  
  igl = 0;        xl = matfunc->xval[igl];
  igr = ngrd - 1; xr = matfunc->xval[igr];
  
  check((xx-xl) * (xx-xr) <= 0.0, "Invalid xx: %15.5e", xx);
  
  while(igr - igl > 1) {
    igm = (igl + igr) / 2; xm = matfunc->xval[igm];
    if((xx-xl) * (xx-xm) <= 0.0) {
      igr = igm; xr = xm;
    } else if((xx-xm) * (xx-xr) <= 0.0) {
      igl = igm; xl = xm;
    }
  }
  
  matl = matfunc->ptr[igl];
  matr = matfunc->ptr[igr];
  for(icol = 0; icol < ncol; ++icol) {
    for(irow = 0; irow < nrow; ++irow) {
      mat->ptr[icol][irow] = matl->ptr[icol][irow]
        + (xx - xl) / (xr - xl) *(matr->ptr[icol][irow] - matl->ptr[icol][irow]);
    }
  }
  
  return;
  
 error:
  abort();
}

void matfunc_intrp_tag_dcmplx(matfunc_dcmplx *matfunc)
{
  int          ngrd, ig, nrow, irow, ncol, icol, ntag, nlen;
  dreal        fre, fim, fremax, fimmax, absmax, dfre, dfim;
  dcmplx       zl, zr;
  mat2d_dcmplx *mat = NULL, *matl = NULL, *matr = NULL;
  
  ngrd = matfunc->ngrd;
  nrow = matfunc->nrow;
  ncol = matfunc->ncol;
  
  check(ngrd > 1, "Invalid ngrd: %5d", ngrd);
  
  nlen = ngrd * sizeof(int);
  memset(matfunc->tag, 0, nlen);
  matfunc->ntag = 0;
  
  fremax = 0.0; fimmax = 0.0;
  for(ig = 0; ig < ngrd; ++ig) {
    mat = matfunc->ptr[ig];
    for(icol = 0; icol < ncol; ++icol) {
      for(irow = 0; irow < nrow; ++irow) {
        fre = fabs(creal(mat->ptr[icol][irow]));
        fim = fabs(cimag(mat->ptr[icol][irow]));
        if(fre > fremax) fremax = fre;
        if(fim > fimmax) fimmax = fim;
      }
    }
  }
  
  absmax = sqrt(fremax*fremax + fimmax*fimmax);
  
  for(ig = 0, ntag = 0; ig < ngrd - 1; ++ig) {
    matl = matfunc->ptr[ig];
    matr = matfunc->ptr[ig+1];
    for(icol = 0; icol < ncol; ++icol) {
      for(irow = 0; irow < nrow; ++irow) {
        zl  = matl->ptr[icol][irow];
        zr  = matr->ptr[icol][irow];
        dfre= creal(zr - zl);
        dfim= cimag(zr - zl);
        if(fabs(dfre) > RATIO * absmax || fabs(dfim) > RATIO * absmax) {
          matfunc->tag[ig] = 1;
        }
      }
    }
    if(matfunc->tag[ig] == 1) ++ntag;
  }
  matfunc->ntag = ntag;
  
  return;
  
 error:
  abort();
}

// For debug
/*
#include <math.h>
#include "lapack.h"

void matfunc_check_dreal(void)
{
  const int     ngrd = 2001;
  const int     ndim = 1;
  const dreal   xmin =-10.0;
  const dreal   xmax = 10.0;
  int           ig, idim, jdim;
  dreal         xx;
  mat2d_dreal   *mat = NULL;
  matfunc_dreal *matfunc = NULL;
  
  matfunc = matfunc_new_dreal();
  check_mem(matfunc, "matfunc");
  matfunc->alloc(ngrd, ndim, ndim, matfunc);
  
  mat = mat2d_new_dreal();
  check_mem(mat, "mat");
  mat->alloc(ndim, ndim, mat);
  
  for(ig = 0; ig < ngrd; ++ig) {
    xx = xmin + ig * (xmax-xmin) / (ngrd - 1);
    matfunc->xval[ig] = xx;
    for(idim = 0; idim < ndim; ++idim) {
      for(jdim = 0; jdim < ndim; ++jdim) {
        matfunc->ptr[ig]->ptr[idim][jdim] = exp(-xx*xx);
      }
    }
  }
  
  xx =-2.618;
  matfunc_bisec_dreal(xx, matfunc, mat);
  printf("%15.5e %15.5e\n", mat->ptr[0][0], exp(-xx*xx));
  
  matfunc_del_dreal(matfunc);
  mat2d_del_dreal(mat);
  
  return;
  
 error:
  abort();
}

void matfunc_check_dcmplx(void)
{
  const int      ngrd = 2001;
  const int      ndim = 1;
  const dreal    xmin =-10.0;
  const dreal    xmax = 10.0;
  int            ig, idim, jdim;
  dreal          xx;
  mat2d_dcmplx   *mat = NULL;
  matfunc_dcmplx *matfunc = NULL;
  
  matfunc = matfunc_new_dcmplx();
  check_mem(matfunc, "matfunc");
  matfunc->alloc(ngrd, ndim, ndim, matfunc);
  
  mat = mat2d_new_dcmplx();
  check_mem(mat, "mat");
  mat->alloc(ndim, ndim, mat);
  
  for(ig = 0; ig < ngrd; ++ig) {
    xx = xmin + ig * (xmax-xmin) / (ngrd - 1);
    matfunc->xval[ig] = xx;
    for(idim = 0; idim < ndim; ++idim) {
      for(jdim = 0; jdim < ndim; ++jdim) {
        matfunc->ptr[ig]->ptr[idim][jdim] = exp(-xx*xx);
      }
    }
  }
  
  xx =-2.618;
  matfunc_bisec_dcmplx(xx, matfunc, mat);
  printf("%15.5e %15.5e\n", creal(mat->ptr[0][0]), exp(-xx*xx));
  
  matfunc_del_dcmplx(matfunc);
  mat2d_del_dcmplx(mat);
  
  return;
  
 error:
  abort();
}
*/
