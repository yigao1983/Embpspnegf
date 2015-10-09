#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "dbg.h"
#include "sys.h"
#include "constants.h"
#include "numeric.h"
#include "variables.h"
#include "lapack.h"
#include "stats.h"
#include "dmft.h"
#include "nca.h"

static const char *ncalabel = "nca:            ";
static const int  MAXITER = 20;

static matfunc_dcmplx **GrnLpsp_aug = NULL, **GrnGpsp_aug = NULL,
                      **GrnRpsp_aug = NULL, **GrnApsp_aug = NULL;

typedef void (*sgmmatfunc)(int, dreal, mat2d_dcmplx *);
typedef void (*grnmatfunc)(int, dreal, const mat2d_dcmplx *, mat2d_dcmplx *);
typedef void (*gqsmatfunc)(int, dreal, mat2d_dcmplx *);

static void nca_out_dospsp(void);
static void nca_out_poppsp(void);

static void nca_sgmrpspmat_init(int iblk, dreal ww, mat2d_dcmplx *SgmRpspmat)
{
  const int    thr = 1e-9;
  int          npspblki, npspblkj, jblk, ipspblki, jpspblki, ipspblkj,
               ipspi, jpspi, ipspj, ispin, iq, jq;
  dreal        wwg, wwl;
  dcmplx       xil, xir;
  dcmplx       lambdag, lambdal;
  mat2d_dcmplx *LamGmat = NULL, *LamLmat = NULL;
  
  SgmRpspmat->reset(SgmRpspmat);
  
  LamGmat = mat2d_new_dcmplx();
  check_mem(LamGmat, "LamGmat");
  LamGmat->alloc(Nqus, Nqus, LamGmat);
  
  LamLmat = mat2d_new_dcmplx();
  check_mem(LamLmat, "LamLmat");
  LamLmat->alloc(Nqus, Nqus, LamLmat);
  
  npspblki = Npspblk[iblk];
  
  for(ipspblki = 0; ipspblki < npspblki; ++ipspblki) {
    ipspi = Pspblk[iblk][ipspblki];
    for(jpspblki = 0; jpspblki < npspblki; ++jpspblki) {
      jpspi = Pspblk[iblk][jpspblki];
      // Find nonvanishing overlap
      for(ispin = 0; ispin < Nspin; ++ispin) {
        for(iq = 0; iq < Nqus; ++iq) {
          for(jq = 0; jq < Nqus; ++jq) {
            
            for(jblk = 0; jblk < Nblk; ++jblk) {
              npspblkj = Npspblk[jblk];
              for(ipspblkj = 0; ipspblkj < npspblkj; ++ipspblkj) {
                ipspj = Pspblk[jblk][ipspblkj];
                
                // Greater Lambda part
                xil = Xipsp[ispin][jq]->ptr[ipspj][jpspi];
                xir = Xipsp[ispin][iq]->ptr[ipspj][ipspi];
                
                if(cabs(xil * conj(xir)) > thr) {
                  wwg = ww - Epsp[ipspj];
                  matfunc_bisec_dcmplx(wwg, LamG_aug[ispin], LamGmat);
                  lambdag = LamGmat->ptr[iq][jq];
                  SgmRpspmat->ptr[ipspblki][jpspblki] += 0.5 * xil * conj(xir) * lambdag;
                }
                
                // Lesser Lambda part
                xil = Xipsp[ispin][iq]->ptr[ipspi][ipspj];
                xir = Xipsp[ispin][jq]->ptr[jpspi][ipspj];
                
                if(cabs(xil * conj(xir)) > thr) {
                  wwl = Epsp[ipspj] - ww;
                  matfunc_bisec_dcmplx(wwl, LamL_aug[ispin], LamLmat);
                  lambdal = LamLmat->ptr[iq][jq];
                  SgmRpspmat->ptr[ipspblki][jpspblki] -= 0.5 * xil * conj(xir) * lambdal;
                }
              }
            }
          }
        }
      }
    }
    SgmRpspmat->ptr[ipspblki][ipspblki] -= I * INFTES;
    if(cimag(SgmRpspmat->ptr[ipspblki][ipspblki]) > 0.0) {
      SgmRpspmat->ptr[ipspblki][ipspblki] =-I * INFTES;
    }
  }
  
  mat2d_del_dcmplx(LamLmat);
  mat2d_del_dcmplx(LamGmat);
  
  return;
  
 error:
  abort();
}

static void nca_sgmrpspmat_iter(int iblk, dreal ww, mat2d_dcmplx *SgmRpspmat)
{
  const int      thr = 1e-9;
  int            ngrd, npspblki, npspblkj, jblk, ipspblki, jpspblki, ipspblkj, jpspblkj,
                 ipspi, jpspi, ipspj, jpspj, ispin, iq, jq, ig;
  dreal          wwg, wwl;
  dcmplx         xil, xir, dSgm;
  dcmplx         lambdag, lambdal;
  size_t         nlen;
  dreal          *xval = NULL;
  dcmplx         *lamxgrn = NULL;
  mat2d_dcmplx   *LamGmat = NULL, *LamLmat = NULL, *GrnRpspmat = NULL;
  matfunc_dcmplx *GrnRpspblk = NULL;
  
  SgmRpspmat->reset(SgmRpspmat);
  
  LamGmat = mat2d_new_dcmplx();
  check_mem(LamGmat, "LamGmat");
  LamGmat->alloc(Nqus, Nqus, LamGmat);
  
  LamLmat = mat2d_new_dcmplx();
  check_mem(LamLmat, "LamLmat");
  LamLmat->alloc(Nqus, Nqus, LamLmat);
  
  npspblki = Npspblk[iblk];
  
  for(ipspblki = 0; ipspblki < npspblki; ++ipspblki) {
    ipspi = Pspblk[iblk][ipspblki];
    for(jpspblki = 0; jpspblki < npspblki; ++jpspblki) {
      jpspi = Pspblk[iblk][jpspblki];
      // Find nonvanishing overlap
      for(ispin = 0; ispin < Nspin; ++ispin) {
        for(iq = 0; iq < Nqus; ++iq) {
          for(jq = 0; jq < Nqus; ++jq) {
            
            for(jblk = 0; jblk < Nblk; ++jblk) {
              npspblkj   = Npspblk[jblk];
              GrnRpspblk = GrnRpsp[jblk];
              
              ngrd    = GrnRpspblk->ngrd;
              xval    = (dreal *) calloc(ngrd, sizeof(dreal));
              lamxgrn = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
              
              nlen = ngrd * sizeof(dreal);
              memcpy(xval, GrnRpspblk->xval, nlen);
              
              for(ipspblkj = 0; ipspblkj < npspblkj; ++ipspblkj) {
                ipspj = Pspblk[jblk][ipspblkj];
                for(jpspblkj = 0; jpspblkj < npspblkj; ++jpspblkj) {
                  jpspj = Pspblk[jblk][jpspblkj];
                  // Greater Lambda part
                  xil = Xipsp[ispin][jq]->ptr[jpspj][jpspi];
                  xir = Xipsp[ispin][iq]->ptr[ipspj][ipspi];
                  
                  if(cabs(xil * conj(xir)) > thr) {
                    for(ig = 0; ig < ngrd; ++ig) {
                      GrnRpspmat = GrnRpspblk->ptr[ig];
                      wwg = ww - xval[ig];
                      matfunc_bisec_dcmplx(wwg, LamG_aug[ispin], LamGmat);
                      lambdag = LamGmat->ptr[iq][jq];
                      lamxgrn[ig] = I / (2.0*PI) * xil * conj(xir)
                        * lambdag * GrnRpspmat->ptr[ipspblkj][jpspblkj];
                    }
                    dSgm = numeric_quadrat_dcmplx(ngrd, xval, lamxgrn);
                    SgmRpspmat->ptr[ipspblki][jpspblki] += dSgm;
                  }
                  
                  // Lesser Lambda part
                  xil = Xipsp[ispin][iq]->ptr[ipspi][ipspj];
                  xir = Xipsp[ispin][jq]->ptr[jpspi][jpspj];
                  
                  if(cabs(xil * conj(xir)) > thr) {
                    for(ig = 0; ig < ngrd; ++ig) {
                      GrnRpspmat = GrnRpspblk->ptr[ig];
                      wwl = xval[ig] - ww;
                      matfunc_bisec_dcmplx(wwl, LamL_aug[ispin], LamLmat);
                      lambdal = LamLmat->ptr[iq][jq];
                      lamxgrn[ig] =-I / (2.0*PI) * xil * conj(xir)
                        * lambdal * GrnRpspmat->ptr[ipspblkj][jpspblkj];
                    }
                    dSgm = numeric_quadrat_dcmplx(ngrd, xval, lamxgrn);
                    SgmRpspmat->ptr[ipspblki][jpspblki] += dSgm;
                  }
                }
              }
              freeup(lamxgrn); freeup(xval);
            }
          }
        }
      }
    }
    SgmRpspmat->ptr[ipspblki][ipspblki] -= I * INFTES;
    if(cimag(SgmRpspmat->ptr[ipspblki][ipspblki]) > 0.0) {
      SgmRpspmat->ptr[ipspblki][ipspblki] =-I * INFTES;
    }
  }
  
  mat2d_del_dcmplx(LamLmat);
  mat2d_del_dcmplx(LamGmat);
  
  return;
  
 error:
  abort();
}

static void nca_sgmlpspmat_init(int iblk, dreal ww, mat2d_dcmplx *SgmLpspmat)
{
  const int    thr = 1e-9;
  int          npspblki, npspblkj, jblk, ipspblki, jpspblki, ipspblkj,
               ipspi, jpspi, ipspj, ispin, iq, jq, zeta;
  dreal        wwl, wwg;
  dcmplx       xil, xir;
  dcmplx       lambdal, lambdag;
  mat2d_dcmplx *LamLmat = NULL, *LamGmat = NULL;
  
  SgmLpspmat->reset(SgmLpspmat);
  
  LamLmat = mat2d_new_dcmplx();
  check_mem(LamLmat, "LamLmat");
  LamLmat->alloc(Nqus, Nqus, LamLmat);
  
  LamGmat = mat2d_new_dcmplx();
  check_mem(LamGmat, "LamGmat");
  LamGmat->alloc(Nqus, Nqus, LamGmat);
  
  npspblki = Npspblk[iblk];
  
  for(ipspblki = 0; ipspblki < npspblki; ++ipspblki) {
    ipspi = Pspblk[iblk][ipspblki];
    for(jpspblki = 0; jpspblki < npspblki; ++jpspblki) {
      jpspi = Pspblk[iblk][jpspblki];
      // Find nonvanishing overlap
      for(ispin = 0; ispin < Nspin; ++ispin) {
        for(iq = 0; iq < Nqus; ++iq) {
          for(jq = 0; jq < Nqus; ++jq) {
            
            for(jblk = 0; jblk < Nblk; ++jblk) {
              npspblkj = Npspblk[jblk];
              for(ipspblkj = 0; ipspblkj < npspblkj; ++ipspblkj) {
                ipspj = Pspblk[jblk][ipspblkj];
                zeta  = Zeta[ipspj];
                // Lesser Lambda part
                xil = Xipsp[ispin][jq]->ptr[ipspj][jpspi];
                xir = Xipsp[ispin][iq]->ptr[ipspj][ipspi];
                
                if(cabs(xil * conj(xir)) > thr) {
                  wwl = ww - Epsp[ipspj];
                  matfunc_bisec_dcmplx(wwl, LamL_aug[ispin], LamLmat);
                  lambdal = LamLmat->ptr[iq][jq];
                  SgmLpspmat->ptr[ipspblki][jpspblki] += zeta * xil * conj(xir) * lambdal / Npsp;
                }
                
                // Greater Lambda part
                xil = Xipsp[ispin][iq]->ptr[ipspi][ipspj];
                xir = Xipsp[ispin][jq]->ptr[jpspi][ipspj];
                
                if(cabs(xil * conj(xir)) > thr) {
                  wwg = Epsp[ipspj] - ww;
                  matfunc_bisec_dcmplx(wwg, LamG_aug[ispin], LamGmat);
                  lambdag = LamGmat->ptr[iq][jq];
                  SgmLpspmat->ptr[ipspblki][jpspblki] -= zeta * xil * conj(xir) * lambdag / Npsp;
                }
              }
            }
          }
        }
      }
    }
    SgmLpspmat->ptr[ipspblki][ipspblki] -= Zeta[ipspi] * I * INFTES;
    if(cabs(SgmLpspmat->ptr[ipspblki][ipspblki]) < INFTES) {
      SgmLpspmat->ptr[ipspblki][ipspblki] =-Zeta[ipspi] * I * INFTES;
    }
  }
  
  mat2d_del_dcmplx(LamGmat);
  mat2d_del_dcmplx(LamLmat);
  
  return;
  
 error:
  abort();
}

static void nca_sgmlpspmat_iter(int iblk, dreal ww, mat2d_dcmplx *SgmLpspmat)
{
  const int      thr = 1e-9;
  int            ngrd, npspblki, npspblkj, jblk, ipspblki, jpspblki, ipspblkj, jpspblkj,
                 ipspi, jpspi, ipspj, jpspj, ispin, iq, jq, ig;
  dreal          wwl, wwg;
  dcmplx         xil, xir, dSgm;
  dcmplx         lambdal, lambdag;
  size_t         nlen;
  dreal          *xval = NULL;
  dcmplx         *lamxgrn = NULL;
  mat2d_dcmplx   *LamLmat = NULL, *LamGmat = NULL, *GrnLpspmat = NULL;
  matfunc_dcmplx *GrnLpspblk = NULL;
  
  SgmLpspmat->reset(SgmLpspmat);
  
  LamLmat = mat2d_new_dcmplx();
  check_mem(LamLmat, "LamLmat");
  LamLmat->alloc(Nqus, Nqus, LamLmat);
  
  LamGmat = mat2d_new_dcmplx();
  check_mem(LamGmat, "LamGmat");
  LamGmat->alloc(Nqus, Nqus, LamGmat);
  
  npspblki = Npspblk[iblk];
  
  for(ipspblki = 0; ipspblki < npspblki; ++ipspblki) {
    ipspi = Pspblk[iblk][ipspblki];
    for(jpspblki = 0; jpspblki < npspblki; ++jpspblki) {
      jpspi = Pspblk[iblk][jpspblki];
      // Find nonvanishing overlap
      for(ispin = 0; ispin < Nspin; ++ispin) {
        for(iq = 0; iq < Nqus; ++iq) {
          for(jq = 0; jq < Nqus; ++jq) {
            
            for(jblk = 0; jblk < Nblk; ++jblk) {
              npspblkj   = Npspblk[jblk];
              GrnLpspblk = GrnLpsp[jblk];
              
              ngrd    = GrnLpspblk->ngrd;
              xval    = (dreal *) calloc(ngrd, sizeof(dreal));
              lamxgrn = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
              
              nlen = ngrd * sizeof(dreal);
              memcpy(xval, GrnLpspblk->xval, nlen);
              
              for(ipspblkj = 0; ipspblkj < npspblkj; ++ipspblkj) {
                ipspj = Pspblk[jblk][ipspblkj];
                for(jpspblkj = 0; jpspblkj < npspblkj; ++jpspblkj) {
                  jpspj = Pspblk[jblk][jpspblkj];
                  // Lesser Lambda part
                  xil = Xipsp[ispin][jq]->ptr[jpspj][jpspi];
                  xir = Xipsp[ispin][iq]->ptr[ipspj][ipspi];
                  
                  if(cabs(xil * conj(xir)) > thr) {
                    for(ig = 0; ig < ngrd; ++ig) {
                      GrnLpspmat = GrnLpspblk->ptr[ig];
                      wwl = ww - xval[ig];
                      matfunc_bisec_dcmplx(wwl, LamL_aug[ispin], LamLmat);
                      lambdal = LamLmat->ptr[iq][jq];
                      lamxgrn[ig] = I / (2.0*PI) * xil * conj(xir) 
                        * lambdal * GrnLpspmat->ptr[ipspblkj][jpspblkj];
                    }
                    dSgm = numeric_quadrat_dcmplx(ngrd, xval, lamxgrn);
                    SgmLpspmat->ptr[ipspblki][jpspblki] += dSgm;
                  }
                  // Greater Lambda part
                  xil = Xipsp[ispin][iq]->ptr[ipspi][ipspj];
                  xir = Xipsp[ispin][jq]->ptr[jpspi][jpspj];
                  
                  if(cabs(xil * conj(xir)) > thr) {
                    for(ig = 0; ig < ngrd; ++ig) {
                      GrnLpspmat = GrnLpspblk->ptr[ig];
                      wwg = xval[ig] - ww;
                      matfunc_bisec_dcmplx(wwg, LamG_aug[ispin], LamGmat);
                      lambdag = LamGmat->ptr[iq][jq];
                      lamxgrn[ig] =-I / (2.0*PI) * xil * conj(xir)
                        * lambdag * GrnLpspmat->ptr[ipspblkj][jpspblkj];
                    }
                    dSgm = numeric_quadrat_dcmplx(ngrd, xval, lamxgrn);
                    SgmLpspmat->ptr[ipspblki][jpspblki] += dSgm;
                  }
                }
              }
              freeup(lamxgrn); freeup(xval);
            }
          }
        }
      }
    }
    SgmLpspmat->ptr[ipspblki][ipspblki] -= Zeta[ipspi] * I * INFTES;
    if(cabs(SgmLpspmat->ptr[ipspblki][ipspblki]) < INFTES) {
      SgmLpspmat->ptr[ipspblki][ipspblki] =-Zeta[ipspi] * I * INFTES;
    }
  }
  
  mat2d_del_dcmplx(LamGmat);
  mat2d_del_dcmplx(LamLmat);
  
  return;
  
 error:
  abort();
}

static void nca_grnrpspmat(int iblk, dreal ww, const mat2d_dcmplx *SgmRpspmat, mat2d_dcmplx *GrnRpspmat)
{
  int    npspblki, ipspblk, jpspblk;
  dcmplx **Hptr = NULL, **Sptr = NULL, **Gptr = NULL;
  
  npspblki = Npspblk[iblk];
  
  Hptr = Hpsp[iblk]->ptr;
  Sptr = SgmRpspmat->ptr;
  Gptr = GrnRpspmat->ptr;
  
  for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
    for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
      Gptr[ipspblk][jpspblk] =-Hptr[ipspblk][jpspblk] - Sptr[ipspblk][jpspblk];
    }
    Gptr[ipspblk][ipspblk] += ww;
  }
  
  lapack_zinvs(npspblki, GrnRpspmat->addr);
}

static void nca_grnlpspmat(int iblk, dreal ww, const mat2d_dcmplx *SgmLpspmat, mat2d_dcmplx *GrnLpspmat)
{
  int            npspblki;
  dcmplx         *mat = NULL;
  mat2d_dcmplx   *GrnRpspmat = NULL;
  matfunc_dcmplx *GrnRpspblk = NULL;
  
  npspblki   = Npspblk[iblk];
  GrnRpspblk = GrnRpsp[iblk];
  
  mat = (dcmplx *) calloc(npspblki*npspblki, sizeof(dcmplx));
  check_mem(mat, "mat");
  
  GrnRpspmat = mat2d_new_dcmplx();
  check_mem(GrnRpspmat, "GrnRpspmat");
  GrnRpspmat->alloc(npspblki, npspblki, GrnRpspmat);
  
  matfunc_bisec_dcmplx(ww, GrnRpspblk, GrnRpspmat);
  
  lapack_zgemm(npspblki, npspblki, npspblki, 'N', 'C',
               1.0, SgmLpspmat->addr, GrnRpspmat->addr,
               0.0, mat);
  lapack_zgemm(npspblki, npspblki, npspblki, 'N', 'N',
               1.0, GrnRpspmat->addr, mat,
               0.0, GrnLpspmat->addr);
  
  mat2d_del_dcmplx(GrnRpspmat);
  freeup(mat);
  
  return;
  
 error:
  abort();
}

static void nca_grnlpsp_renorm(void)
{
  int            ipsp, iblk, ngrd, ig, npspblki, ipspblk, jpspblk;
  dreal          totpop;
  mat2d_dcmplx   *GrnLpspmat = NULL;
  matfunc_dcmplx *GrnLpspblk = NULL;
  
  for(ipsp = 0, totpop = 0.0; ipsp < Npsp; ++ipsp) {
    totpop += Sumpoppsp[ipsp];
  }
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki   = Npspblk[iblk];
    GrnLpspblk = GrnLpsp[iblk];
    ngrd = GrnLpspblk->ngrd;
    for(ig = 0; ig < ngrd; ++ig) {
      GrnLpspmat = GrnLpspblk->ptr[ig];
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
          GrnLpspmat->ptr[ipspblk][jpspblk] /= totpop;
        }
      }
    }
  }
}

static void nca_dospsp(void)
{
  int            iblk, npspblki, ipspblk, ipsp, ngrd, ig;
  size_t         nlen;
  dreal          *xval = NULL, *dospspi = NULL;
  mat2d_dcmplx   *GrnRpspmat = NULL;
  matfunc_dcmplx *GrnRpspblk = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki = Npspblk[iblk];
    GrnRpspblk = GrnRpsp[iblk];
    ngrd = GrnRpspblk->ngrd;
    nlen = ngrd * sizeof(dreal);
    xval = (dreal *) calloc(ngrd, sizeof(dreal));
    check_mem(xval, "xval");
    memcpy(xval, GrnRpspblk->xval, nlen);
    for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
      ipsp = Pspblk[iblk][ipspblk];
      dospspi = Dospsp[ipsp];
      if(dospspi) freeup(Dospsp[ipsp]);
      dospspi = (dreal *) calloc(ngrd, sizeof(dreal));
      check_mem(dospspi, "dospspi");
      for(ig = 0; ig < ngrd; ++ig) {
        GrnRpspmat = GrnRpspblk->ptr[ig];
        dospspi[ig] =-cimag(GrnRpspmat->ptr[ipspblk][ipspblk]) / PI;
      }
      Sumdospsp[ipsp] = numeric_quadrat_dreal(ngrd, xval, dospspi);
      Dospsp[ipsp] = dospspi;
    }
    freeup(xval);
  }
  
  return;
  
 error:
  if(xval)    freeup(xval);
  if(dospspi) freeup(dospspi);
  abort();
}

static void nca_poppsp(void)
{
  int            iblk, npspblki, ipspblk, ipsp, ngrd, ig, zeta;
  size_t         nlen;
  dreal          *xval = NULL, *poppspi = NULL;
  mat2d_dcmplx   *GrnLpspmat = NULL;
  matfunc_dcmplx *GrnLpspblk = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki   = Npspblk[iblk];
    GrnLpspblk = GrnLpsp[iblk];
    ngrd = GrnLpspblk->ngrd;
    nlen = ngrd * sizeof(dreal);
    xval = (dreal *) calloc(ngrd, sizeof(dreal));
    check_mem(xval, "xval");
    memcpy(xval, GrnLpspblk->xval, nlen);
    for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
      ipsp = Pspblk[iblk][ipspblk];
      zeta = Zeta[ipsp];
      poppspi = Poppsp[ipsp];
      if(poppspi) freeup(Poppsp[ipsp]);
      poppspi = (dreal *) calloc(ngrd, sizeof(dreal));
      check_mem(poppspi, "poppspi");
      for(ig = 0; ig < ngrd; ++ig) {
        GrnLpspmat = GrnLpspblk->ptr[ig];
        poppspi[ig] =-zeta * cimag(GrnLpspmat->ptr[ipspblk][ipspblk]) / (2.0*PI);
      }
      Sumpoppsp[ipsp] = numeric_quadrat_dreal(ngrd, xval, poppspi);
      Poppsp[ipsp] = poppspi;
    }
    freeup(xval);
  }
  
  return;
  
 error:
  if(xval)    freeup(xval);
  if(poppspi) freeup(poppspi);
  abort();
}

static void nca_psp_intrp(int iblk, const matfunc_dcmplx *Sgmpspblk, const matfunc_dcmplx *Grnpspblk,
                          matfunc_dcmplx *Sgmpspblk_aug, matfunc_dcmplx *Grnpspblk_aug,
                          sgmmatfunc nca_sgmpspmat, grnmatfunc nca_grnpspmat)
{
  int          ngrd, ig, ngrd_aug, ig_aug, ntag, nrow, ncol, isub;
  size_t       nlen;
  dreal        dxval, dxsub, xsub, xsubmin;
  dcmplx       *Sgmbuf = NULL, *Grnbuf = NULL;
  mat2d_dcmplx *Grnpspmat = NULL, *Grnpspmat_aug = NULL,
               *Sgmpspmat = NULL, *Sgmpspmat_aug = NULL;
  
  ngrd = Grnpspblk->ngrd;
  ntag = Grnpspblk->ntag;
  ngrd_aug = ngrd + (NSUB - 1) * ntag;
  
  nrow = Grnpspblk->nrow;
  ncol = Grnpspblk->ncol;
  
  nlen = sizeof(dcmplx)*nrow*ncol;
  
  Sgmpspblk_aug->deall(Sgmpspblk_aug);
  Grnpspblk_aug->deall(Grnpspblk_aug);
  
  Sgmpspblk_aug->alloc(ngrd_aug, nrow, ncol, Sgmpspblk_aug);
  Grnpspblk_aug->alloc(ngrd_aug, nrow, ncol, Grnpspblk_aug);
  
  Sgmbuf = (dcmplx *) calloc(nrow*ncol*ngrd_aug, sizeof(dcmplx));
  Grnbuf = (dcmplx *) calloc(nrow*ncol*ngrd_aug, sizeof(dcmplx));
  
  for(ig = 0, ig_aug = 0; ig < ngrd - 1; ++ig, ++ig_aug) {
    Grnpspblk_aug->xval[ig_aug] = Grnpspblk->xval[ig];
    Sgmpspblk_aug->xval[ig_aug] = Sgmpspblk->xval[ig];
    if(Grnpspblk->tag[ig]) {
      dxval = Grnpspblk->xval[ig+1] - Grnpspblk->xval[ig];
      dxsub = dxval / NSUB;
      xsubmin = Grnpspblk->xval[ig];
      for(isub = 1; isub < NSUB; ++isub) {
        ++ig_aug;
        xsub = xsubmin + dxsub * isub;
        Grnpspblk_aug->xval[ig_aug] = xsub;
        Sgmpspblk_aug->xval[ig_aug] = xsub;
      }
    }
  }
  
  for(ig = 0, ig_aug = 0; ig < ngrd - 1; ++ig, ++ig_aug) {
    Grnpspmat = Grnpspblk->ptr[ig];
    Sgmpspmat = Sgmpspblk->ptr[ig];
    Grnpspmat_aug = Grnpspblk_aug->ptr[ig_aug];
    Sgmpspmat_aug = Sgmpspblk_aug->ptr[ig_aug];
    
    if(ig_aug % Size == Rank) {
      mat2d_copy_dcmplx(Grnpspmat, Grnpspmat_aug);
      mat2d_copy_dcmplx(Sgmpspmat, Sgmpspmat_aug);
      memcpy(Grnbuf + nrow*ncol*ig_aug, Grnpspmat_aug->addr, nlen);
      memcpy(Sgmbuf + nrow*ncol*ig_aug, Sgmpspmat_aug->addr, nlen);
    }
    if(Grnpspblk->tag[ig]) {
      for(isub = 1; isub < NSUB; ++isub) {
        ++ig_aug;
        xsub = Grnpspblk_aug->xval[ig_aug];
        Grnpspmat_aug = Grnpspblk_aug->ptr[ig_aug];
        Sgmpspmat_aug = Sgmpspblk_aug->ptr[ig_aug];
        if(ig_aug % Size == Rank) {
          nca_sgmpspmat(iblk, xsub, Sgmpspmat_aug);
          nca_grnpspmat(iblk, xsub, Sgmpspmat_aug, Grnpspmat_aug);
          memcpy(Grnbuf + nrow*ncol*ig_aug, Grnpspmat_aug->addr, nlen);
          memcpy(Sgmbuf + nrow*ncol*ig_aug, Sgmpspmat_aug->addr, nlen);
        }
      }
    }
  }
  
  sys_allreduce_sum_dcmplx(Sgmbuf, nrow*ncol*ngrd_aug); sys_sync();
  sys_allreduce_sum_dcmplx(Grnbuf, nrow*ncol*ngrd_aug); sys_sync();
  
  for(ig_aug = 0; ig_aug < ngrd_aug - 1; ++ig_aug) {
    Sgmpspmat_aug = Sgmpspblk_aug->ptr[ig_aug];
    Grnpspmat_aug = Grnpspblk_aug->ptr[ig_aug];
    if(ig_aug % Size != Rank) {
      memcpy(Sgmpspmat_aug->addr, Sgmbuf + nrow*ncol*ig_aug, nlen);
      memcpy(Grnpspmat_aug->addr, Grnbuf + nrow*ncol*ig_aug, nlen);
    }
  }
  
  Grnpspblk_aug->xval[ngrd_aug-1] = Grnpspblk->xval[ngrd-1];
  mat2d_copy_dcmplx(Grnpspblk->ptr[ngrd-1], Grnpspblk_aug->ptr[ngrd_aug-1]);
  Sgmpspblk_aug->xval[ngrd_aug-1] = Sgmpspblk->xval[ngrd-1];
  mat2d_copy_dcmplx(Sgmpspblk->ptr[ngrd-1], Sgmpspblk_aug->ptr[ngrd_aug-1]);
  
  freeup(Grnbuf); freeup(Sgmbuf);
}

static void nca_psp_unifm(int iblk, matfunc_dcmplx *Sgmpspblk, matfunc_dcmplx *Grnpspblk,
                          sgmmatfunc nca_sgmpspmat, grnmatfunc nca_grnpspmat)
{
  int          npspblki, ig; 
  size_t       nlen;
  dreal        Eng, dE;
  dcmplx       *Sgmbuf = NULL, *Grnbuf = NULL;
  mat2d_dcmplx *Sgmpspmat = NULL, *Grnpspmat = NULL;
  
  npspblki = Npspblk[iblk];
  
  nlen = sizeof(dcmplx)*npspblki*npspblki;
  
  Sgmpspblk->deall(Sgmpspblk);
  Grnpspblk->deall(Grnpspblk);
  
  Sgmpspblk->alloc(Nmesh, npspblki, npspblki, Sgmpspblk);
  Grnpspblk->alloc(Nmesh, npspblki, npspblki, Grnpspblk);
  
  Sgmbuf = (dcmplx *) calloc(npspblki*npspblki*Nmesh, sizeof(dcmplx));
  Grnbuf = (dcmplx *) calloc(npspblki*npspblki*Nmesh, sizeof(dcmplx));
  
  dE = (Emax - Emin) / (Nmesh - 1);
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    Sgmpspblk->xval[ig] = Eng;
    Grnpspblk->xval[ig] = Eng;
  }
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    Sgmpspmat = Sgmpspblk->ptr[ig];
    Grnpspmat = Grnpspblk->ptr[ig];
    
    if(ig % Size == Rank) {
      nca_sgmpspmat(iblk, Eng, Sgmpspmat);
      nca_grnpspmat(iblk, Eng, Sgmpspmat, Grnpspmat);
      memcpy(Sgmbuf + npspblki*npspblki*ig, Sgmpspmat->addr, nlen);
      memcpy(Grnbuf + npspblki*npspblki*ig, Grnpspmat->addr, nlen);
    }
  }
  
  sys_allreduce_sum_dcmplx(Sgmbuf, npspblki*npspblki*Nmesh); sys_sync();
  sys_allreduce_sum_dcmplx(Grnbuf, npspblki*npspblki*Nmesh); sys_sync();
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Sgmpspmat = Sgmpspblk->ptr[ig];
    Grnpspmat = Grnpspblk->ptr[ig];
    if(ig % Size != Rank) {
      memcpy(Sgmpspmat->addr, Sgmbuf + npspblki*npspblki*ig, nlen);
      memcpy(Grnpspmat->addr, Grnbuf + npspblki*npspblki*ig, nlen);
    }
  }
  
  freeup(Grnbuf); freeup(Sgmbuf);
}

static void nca_psp_radpt(sgmmatfunc nca_sgmrpspmat, grnmatfunc nca_grnrpspmat)
{
  int            iblk, cnvg;
  matfunc_dcmplx *SgmRpspblk = NULL, *GrnRpspblk = NULL,
                 *SgmRpspblk_aug = NULL, *GrnRpspblk_aug = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    SgmRpspblk = matfunc_new_dcmplx(); check_mem(SgmRpspblk, "SgmRpspblk");
    GrnRpspblk = matfunc_new_dcmplx(); check_mem(GrnRpspblk, "GrnRpspblk");
    nca_psp_unifm(iblk, SgmRpspblk, GrnRpspblk, nca_sgmrpspmat, nca_grnrpspmat);
    // Adaptive mesh loop
    for(cnvg = 0; GrnRpspblk->ngrd < MAXNGRD && !cnvg; ) {
      //nca_intrp_tag(GrnRpspblk);
      matfunc_intrp_tag_dcmplx(GrnRpspblk);
      if(GrnRpspblk->ntag) {
        SgmRpspblk_aug = matfunc_new_dcmplx(); check_mem(SgmRpspblk_aug, "SgmRpspblk_aug");
        GrnRpspblk_aug = matfunc_new_dcmplx(); check_mem(GrnRpspblk_aug, "GrnRpspblk_aug");
        nca_psp_intrp(iblk, SgmRpspblk, GrnRpspblk, SgmRpspblk_aug, GrnRpspblk_aug,
                      nca_sgmrpspmat, nca_grnrpspmat);
        matfunc_del_dcmplx(SgmRpspblk); SgmRpspblk = SgmRpspblk_aug;
        matfunc_del_dcmplx(GrnRpspblk); GrnRpspblk = GrnRpspblk_aug;
      } else {
        cnvg = 1;
      }
    }
    
    if(Rank == Root) {
      printf("%s iblk, ngrd (R) = %5d %5d\n", ncalabel, iblk, GrnRpspblk->ngrd);
      fflush(stdout);
    }
    matfunc_del_dcmplx(SgmRpsp[iblk]); SgmRpsp[iblk] = SgmRpspblk;
    matfunc_del_dcmplx(GrnRpsp[iblk]); GrnRpsp[iblk] = GrnRpspblk;
  }
  
  return;
  
 error:
  abort();
}

static void nca_psp_ladpt(sgmmatfunc nca_sgmlpspmat, grnmatfunc nca_grnlpspmat)
{
  int            iblk, cnvg;
  matfunc_dcmplx *SgmLpspblk = NULL, *GrnLpspblk = NULL,
                 *SgmLpspblk_aug = NULL, *GrnLpspblk_aug = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    SgmLpspblk = matfunc_new_dcmplx(); check_mem(SgmLpspblk, "SgmLpspblk");
    GrnLpspblk = matfunc_new_dcmplx(); check_mem(GrnLpspblk, "GrnLpspblk");
    nca_psp_unifm(iblk, SgmLpspblk, GrnLpspblk, nca_sgmlpspmat, nca_grnlpspmat);
    // Adaptive mesh loop
    for(cnvg = 0; GrnLpspblk->ngrd < MAXNGRD && !cnvg; ) {
      //nca_intrp_tag(GrnLpspblk);
      matfunc_intrp_tag_dcmplx(GrnLpspblk);
      if(GrnLpspblk->ntag) {
        SgmLpspblk_aug = matfunc_new_dcmplx(); check_mem(SgmLpspblk_aug, "SgmLpspblk_aug");
        GrnLpspblk_aug = matfunc_new_dcmplx(); check_mem(GrnLpspblk_aug, "GrnLpspblk_aug");
        nca_psp_intrp(iblk, SgmLpspblk, GrnLpspblk, SgmLpspblk_aug, GrnLpspblk_aug,
                      nca_sgmlpspmat, nca_grnlpspmat);
        matfunc_del_dcmplx(SgmLpspblk); SgmLpspblk = SgmLpspblk_aug;
        matfunc_del_dcmplx(GrnLpspblk); GrnLpspblk = GrnLpspblk_aug;
      } else {
        cnvg = 1;
      }
    }
    if(Rank == Root) {
      printf("%s iblk, ngrd (L) = %5d %5d\n", ncalabel, iblk, GrnLpspblk->ngrd);
      fflush(stdout);
    }
    matfunc_del_dcmplx(SgmLpsp[iblk]); SgmLpsp[iblk] = SgmLpspblk;
    matfunc_del_dcmplx(GrnLpsp[iblk]); GrnLpsp[iblk] = GrnLpspblk;
  }
  
  return;
  
 error:
  abort();
}

void nca_psp_rscf(void)
{
  int   iter, ipsp, cnvg, nlen;
  dreal dsumdos;
  dreal *sumdos_old = NULL, *sumdos_new = NULL;
  
  if(Rank == Root) {
    printf("%s\n", ncalabel);
    printf("%s ========= SCF Retarded =========\n", ncalabel);
    fflush(stdout);
  }
  
  sumdos_old = (dreal *) calloc(Npsp, sizeof(dreal));
  check_mem(sumdos_old, "sumdos_old");
  sumdos_new = (dreal *) calloc(Npsp, sizeof(dreal));
  check_mem(sumdos_new, "sumdos_new");
  
  nlen = Npsp * sizeof(dreal);
  
  nca_psp_radpt(nca_sgmrpspmat_init, nca_grnrpspmat); nca_dospsp();
  
  for(iter = 0, cnvg = 0; iter < MAXITER && !cnvg; ++iter) {
    
    if(Rank == Root) {
      printf("%s Iter %5d\n", ncalabel, iter); 
      fflush(stdout);
    }
    
    memcpy(sumdos_old, Sumdospsp, nlen);
    
    nca_psp_radpt(nca_sgmrpspmat_iter, nca_grnrpspmat); nca_dospsp();
    
    memcpy(sumdos_new, Sumdospsp, nlen);
    
    if(Rank == Root) {
      printf("%s ipsp, sumdos, dsumdos =\n", ncalabel);
      fflush(stdout);
    }
    for(ipsp = 0, cnvg = 1; ipsp < Npsp; ++ipsp) {
      dsumdos = fabs(sumdos_old[ipsp]-sumdos_new[ipsp]);
      if(dsumdos > TOLDIFF) cnvg = 0; 
      if(Rank == Root) {
        printf("%s %5d %15.5e %15.5e\n", ncalabel, ipsp, sumdos_new[ipsp], dsumdos);
        fflush(stdout);
      }
    }
    
    if(Rank == Root) {
      printf("%s --------------------------------\n", ncalabel);
      fflush(stdout);
    }
  }
  
  freeup(sumdos_new);
  freeup(sumdos_old);
  
  return;
  
 error:
  if(sumdos_old) freeup(sumdos_old);
  if(sumdos_new) freeup(sumdos_new);
  abort();
}

void nca_psp_lscf(void)
{
  int         iter, ipsp, cnvg, nlen;
  dreal       dsumpop, totpop;
  dreal       *sumpop_old = NULL, *sumpop_new = NULL;
  
  if(Rank == Root) {
    printf("%s\n", ncalabel);
    printf("%s ========== SCF Lesser ==========\n", ncalabel);
    fflush(stdout);
  }
  
  sumpop_old = (dreal *) calloc(Npsp, sizeof(dreal));
  check_mem(sumpop_old, "sumpop_old");
  sumpop_new = (dreal *) calloc(Npsp, sizeof(dreal));
  check_mem(sumpop_new, "sumpop_new");
  
  nlen = Npsp * sizeof(dreal);
  
  nca_psp_ladpt(nca_sgmlpspmat_init, nca_grnlpspmat);
  nca_poppsp();
  nca_grnlpsp_renorm();
  
  for(iter = 0, cnvg = 0; iter < MAXITER && !cnvg; ++iter) {
    
    if(Rank == Root) {
      printf("%s Iter %5d\n", ncalabel, iter); 
      fflush(stdout);
    }
    
    memcpy(sumpop_old, Sumpoppsp, nlen);
    
    nca_psp_ladpt(nca_sgmlpspmat_iter, nca_grnlpspmat); nca_poppsp();
    nca_grnlpsp_renorm();
    
    memcpy(sumpop_new, Sumpoppsp, nlen);
    
    if(Rank == Root) {
      printf("%s ipsp, sumpop, dsumpop =\n", ncalabel);
      fflush(stdout);
    }
    for(ipsp = 0, cnvg = 1; ipsp < Npsp; ++ipsp) {
      dsumpop = fabs(sumpop_old[ipsp]-sumpop_new[ipsp]);
      if(dsumpop > TOLDIFF) cnvg = 0; 
      if(Rank == Root) {
        printf("%s %5d %15.5e %15.5e\n", ncalabel, ipsp, sumpop_new[ipsp], dsumpop);
        fflush(stdout);
      }
    }
    for(ipsp = 0, totpop = 0.0; ipsp < Npsp; ++ipsp) {
      totpop += Sumpoppsp[ipsp];
    }
    if(Rank == Root) {
      printf("%s   tot:%15.5e\n", ncalabel, totpop);
      fflush(stdout);
    }
    if(Rank == Root) {
      printf("%s --------------------------------\n", ncalabel);
      fflush(stdout);
    }
  }
  
  freeup(sumpop_new);
  freeup(sumpop_old);
  
  return;
  
 error:
  if(sumpop_old) freeup(sumpop_old);
  if(sumpop_new) freeup(sumpop_new);
  abort();
}

void nca_psp_aprj(void)
{
  int            iblk, npspblki, ngrd, ig, ipspblk, jpspblk;
  size_t         nlen;
  mat2d_dcmplx   *GrnRpspmat = NULL, *SgmRpspmat = NULL,
                 *GrnApspmat = NULL, *SgmApspmat = NULL;
  matfunc_dcmplx *GrnRpspblk = NULL, *SgmRpspblk = NULL,
                 *GrnApspblk = NULL, *SgmApspblk = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki   = Npspblk[iblk];
    GrnRpspblk = GrnRpsp[iblk];
    SgmRpspblk = SgmRpsp[iblk];
    
    GrnApspblk = matfunc_new_dcmplx();
    check_mem(GrnApspblk, "GrnApspblk");
    ngrd = GrnRpspblk->ngrd;
    GrnApspblk->alloc(ngrd, npspblki, npspblki, GrnApspblk);
    nlen = ngrd * sizeof(dreal);
    memcpy(GrnApspblk->xval, GrnRpspblk->xval, nlen);
    
    for(ig = 0; ig < ngrd; ++ig) {
      GrnRpspmat = GrnRpspblk->ptr[ig];
      GrnApspmat = GrnApspblk->ptr[ig];
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
          GrnApspmat->ptr[ipspblk][jpspblk] = conj(GrnRpspmat->ptr[jpspblk][ipspblk]);
        }
      }
    }
    matfunc_del_dcmplx(GrnApsp[iblk]); GrnApsp[iblk] = GrnApspblk;
    
    SgmApspblk = matfunc_new_dcmplx();
    check_mem(SgmApspblk, "SgmApspblk");
    ngrd = SgmRpspblk->ngrd;
    SgmApspblk->alloc(ngrd, npspblki, npspblki, SgmApspblk);
    nlen = ngrd * sizeof(dreal);
    memcpy(SgmApspblk->xval, SgmRpspblk->xval, nlen);
    
    for(ig = 0; ig < ngrd; ++ig) {
      SgmRpspmat = SgmRpspblk->ptr[ig];
      SgmApspmat = SgmApspblk->ptr[ig];
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
          SgmApspmat->ptr[ipspblk][jpspblk] = conj(SgmRpspmat->ptr[jpspblk][ipspblk]);
        }
      }
    }
    matfunc_del_dcmplx(SgmApsp[iblk]); SgmApsp[iblk] = SgmApspblk;
  }
  
  return;
  
 error:
  abort();
}

void nca_psp_gprj(void)
{
  int            iblk, npspblki, ngrd, ig, ipspblk, jpspblk;
  size_t         nlen;
  mat2d_dcmplx   *GrnRpspmat = NULL, *SgmRpspmat = NULL,
                 *GrnGpspmat = NULL, *SgmGpspmat = NULL;
  matfunc_dcmplx *GrnRpspblk = NULL, *SgmRpspblk = NULL,
                 *GrnGpspblk = NULL, *SgmGpspblk = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki   = Npspblk[iblk];
    
    GrnRpspblk = GrnRpsp[iblk];
    SgmRpspblk = SgmRpsp[iblk];
    
    GrnGpspblk = matfunc_new_dcmplx();
    check_mem(GrnGpspblk, "GrnGpspblk");
    ngrd = GrnRpspblk->ngrd;
    GrnGpspblk->alloc(ngrd, npspblki, npspblki, GrnGpspblk);
    nlen = ngrd * sizeof(dreal);
    memcpy(GrnGpspblk->xval, GrnRpspblk->xval, nlen);
    
    for(ig = 0; ig < ngrd; ++ig) {
      GrnRpspmat = GrnRpspblk->ptr[ig];
      GrnGpspmat = GrnGpspblk->ptr[ig];
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
          GrnGpspmat->ptr[ipspblk][jpspblk] = GrnRpspmat->ptr[ipspblk][jpspblk]
            - conj(GrnRpspmat->ptr[jpspblk][ipspblk]);
        }
      }
    }
    matfunc_del_dcmplx(GrnGpsp[iblk]); GrnGpsp[iblk] = GrnGpspblk;
    
    SgmGpspblk = matfunc_new_dcmplx();
    check_mem(SgmGpspblk, "SgmGpspblk");
    ngrd = SgmRpspblk->ngrd;
    SgmGpspblk->alloc(ngrd, npspblki, npspblki, SgmGpspblk);
    nlen = ngrd * sizeof(dreal);
    memcpy(SgmGpspblk->xval, SgmRpspblk->xval, nlen);
    
    for(ig = 0; ig < ngrd; ++ig) {
      SgmRpspmat = SgmRpspblk->ptr[ig];
      SgmGpspmat = SgmGpspblk->ptr[ig];
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
          SgmGpspmat->ptr[ipspblk][jpspblk] = SgmRpspmat->ptr[ipspblk][jpspblk]
            - conj(SgmRpspmat->ptr[jpspblk][ipspblk]);
        }
      }
    }
    matfunc_del_dcmplx(SgmGpsp[iblk]); SgmGpsp[iblk] = SgmGpspblk;
  }
  
  return;
  
 error:
  abort();
}

static void nca_out_dospsp(void)
{
  const char     *dospsp = "Dospsp";
  FILE           *fp = NULL;
  int            iblk, npspblki, ipspblk, ipsp, ngrd, ig;
  char           filename[SHRT_MAX];
  dreal          *xval = NULL;
  matfunc_dcmplx *GrnRpspblk = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki   = Npspblk[iblk];
    GrnRpspblk = GrnRpsp[iblk];
    
    ngrd = GrnRpspblk->ngrd;
    xval = (dreal *) calloc(ngrd, sizeof(dreal));
    check_mem(xval, "xval");
    memcpy(xval, GrnRpspblk->xval, ngrd*sizeof(dreal));
    for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
      ipsp = Pspblk[iblk][ipspblk];
      memset(filename, 0, SHRT_MAX*sizeof(char));
      sprintf(filename, "%s_%04d.ascii", dospsp, ipsp);
      
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", xval[ig], Dospsp[ipsp][ig]);
        }
        fclose(fp);
      }
    }
    freeup(xval);
  }
  
  return;
  
 error:
  abort();
}

static void nca_out_poppsp(void)
{
  const char     *poppsp = "Poppsp";
  FILE           *fp = NULL;
  int            iblk, npspblki, ipspblk, ipsp, ngrd, ig;
  char           filename[SHRT_MAX];
  dreal          *xval = NULL;
  matfunc_dcmplx *GrnLpspblk = NULL;
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki   = Npspblk[iblk];
    GrnLpspblk = GrnLpsp[iblk];
    
    ngrd = GrnLpspblk->ngrd;
    xval = (dreal *) calloc(ngrd, sizeof(dreal));
    check_mem(xval, "xval");
    memcpy(xval, GrnLpspblk->xval, ngrd*sizeof(dreal));
    for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
      ipsp = Pspblk[iblk][ipspblk];
      memset(filename, 0, SHRT_MAX*sizeof(char));
      sprintf(filename, "%s_%04d.ascii", poppsp, ipsp);
      
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", xval[ig], Poppsp[ipsp][ig]);
        }
        fclose(fp);
      }
    }
    freeup(xval);
  }
  
  return;
  
 error:
  abort();
}

void nca_psp_out(void)
{
  nca_out_dospsp();
  nca_out_poppsp();
}

// Following are post-scf calculations for quasiparticle Green's functions
static void nca_matfunc_aug(matfunc_dcmplx *matfunc, matfunc_dcmplx *matfunc_aug)
{
  int   nrow, ncol, ngrd, ngrd_aug, ngrd_neg, ngrd_pos, ig_aug, ig;
  dreal *xval = NULL;
  
  ngrd = matfunc->ngrd;
  ngrd_aug = NAUG * (ngrd - 1) + 1;
  
  nrow = matfunc->nrow;
  ncol = matfunc->ncol;
  
  matfunc_aug->deall(matfunc_aug);
  matfunc_aug->alloc(ngrd_aug, nrow, ncol, matfunc_aug);
  
  xval = matfunc->xval;
  
  for(ig = 0, ngrd_neg = 0, ngrd_pos = 0; ig < ngrd; ++ig) {
    if(xval[ig] > Emin && xval[ig] <= 0.0) ++ngrd_neg;
    if(xval[ig] >= 0.0 && xval[ig] < Emax) ++ngrd_pos;
  }
  
  check(ngrd_neg+ngrd_pos == ngrd-1,
        "Division failure: ngrd_neg+ngrd_pos, ngrd - 1 = %5d %5d\n",
        ngrd_neg+ngrd_pos, ngrd - 1);
  
  for(ig_aug = 0; ig_aug < ngrd_pos; ++ig_aug) {
    ig = ig_aug + ngrd_neg;
    matfunc_aug->xval[ig_aug] = matfunc->xval[ig] + Emin_aug;
    mat2d_copy_dcmplx(matfunc->ptr[ig], matfunc_aug->ptr[ig_aug]);
  }
  
  for(ig_aug = ngrd_pos; ig_aug < ngrd_pos+ngrd; ++ig_aug) {
    ig = ig_aug - ngrd_pos;
    matfunc_aug->xval[ig_aug] = matfunc->xval[ig];
    mat2d_copy_dcmplx(matfunc->ptr[ig], matfunc_aug->ptr[ig_aug]);
  }
  
  for(ig_aug = ngrd_pos+ngrd; ig_aug < ngrd_aug; ++ig_aug) {
    ig = ig_aug - ngrd_pos - ngrd + 1;
    matfunc_aug->xval[ig_aug] = matfunc->xval[ig] + (Emax - Emin);
    mat2d_copy_dcmplx(matfunc->ptr[ig], matfunc_aug->ptr[ig_aug]);
  }
  
  return;
  
 error:
  abort();
}

static void nca_psp_aug_alloc(void)
{
  int iblk;
  
  GrnLpsp_aug = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    GrnLpsp_aug[iblk] = matfunc_new_dcmplx();
    check_mem(GrnLpsp_aug[iblk], "GrnLpsp_aug[iblk]");
    //nca_matfunc_aug(GrnLpsp[iblk], GrnLpsp_aug[iblk]);
    dmft_matfunc_aug(GrnLpsp[iblk], GrnLpsp_aug[iblk]);
  }
  
  GrnGpsp_aug = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    GrnGpsp_aug[iblk] = matfunc_new_dcmplx();
    check_mem(GrnGpsp_aug[iblk], "GrnGpsp_aug[iblk]");
    nca_matfunc_aug(GrnGpsp[iblk], GrnGpsp_aug[iblk]);
  }
  
  GrnRpsp_aug = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    GrnRpsp_aug[iblk] = matfunc_new_dcmplx();
    check_mem(GrnRpsp_aug[iblk], "GrnRpsp_aug[iblk]");
    nca_matfunc_aug(GrnRpsp[iblk], GrnRpsp_aug[iblk]);
  }
  
  GrnApsp_aug = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    GrnApsp_aug[iblk] = matfunc_new_dcmplx();
    check_mem(GrnApsp_aug[iblk], "GrnApsp_aug[iblk]");
    nca_matfunc_aug(GrnApsp[iblk], GrnApsp_aug[iblk]);
  }
  
  return;
  
 error:
  abort();
}

static void nca_psp_aug_deall(void)
{
  int iblk;
  
  if(GrnApsp_aug) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      if(GrnApsp_aug[iblk]) matfunc_del_dcmplx(GrnApsp_aug[iblk]);
    }
    freeup(GrnApsp_aug);
  }
  if(GrnRpsp_aug) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      if(GrnRpsp_aug[iblk]) matfunc_del_dcmplx(GrnRpsp_aug[iblk]);
    }
    freeup(GrnRpsp_aug);
  }
  if(GrnGpsp_aug) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      if(GrnGpsp_aug[iblk]) matfunc_del_dcmplx(GrnGpsp_aug[iblk]);
    }
    freeup(GrnGpsp_aug);
  }
  if(GrnLpsp_aug) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      if(GrnLpsp_aug[iblk]) matfunc_del_dcmplx(GrnLpsp_aug[iblk]);
    }
    freeup(GrnLpsp_aug);
  }
}

static void nca_gxgmat(matfunc_dcmplx *Gfuncj, matfunc_dcmplx *Gfunci,
                       dreal ww, dreal xx, mat2d_dcmplx *GxGmat)
{
  mat2d_dcmplx *Gmatj = NULL, *Gmati = NULL;
  
  Gmatj = mat2d_new_dcmplx();
  check_mem(Gmatj, "Gmatj");
  Gmatj->alloc(NVIB, NVIB, Gmatj);
  
  Gmati = mat2d_new_dcmplx();
  check_mem(Gmati, "Gmati");
  Gmati->alloc(NVIB, NVIB, Gmati);
  
  matfunc_bisec_dcmplx(ww + xx, Gfuncj, Gmatj);
  matfunc_bisec_dcmplx(xx, Gfunci, Gmati);
  lapack_zgemm(NVIB, NVIB, NVIB, 'N', 'N',
               1.0, Gmatj->addr, Gmati->addr,
               0.0, GxGmat->addr);
  
  mat2d_del_dcmplx(Gmati);
  mat2d_del_dcmplx(Gmatj);
  
  return;
  
 error:
  abort();
}

static void nca_gxgfunc_unifm(matfunc_dcmplx *Gfuncj, matfunc_dcmplx *Gfunci,
                              dreal ww, matfunc_dcmplx *GxGfunc)
{
  int          ig;
  dreal        dE, xx;
  mat2d_dcmplx *GxGmat = NULL;
  
  GxGfunc->deall(GxGfunc);
  GxGfunc->alloc(Nmesh, NVIB, NVIB, GxGfunc);
  
  dE = (Emax - Emin) / (Nmesh - 1);
  for(ig = 0; ig < Nmesh; ++ig) {
    xx = Emin + dE * ig;
    GxGfunc->xval[ig] = xx;
    GxGmat = GxGfunc->ptr[ig];
    nca_gxgmat(Gfuncj, Gfunci, ww, xx, GxGmat);
  }
}

static void nca_gxgfunc_intrp(matfunc_dcmplx *Gfuncj, matfunc_dcmplx *Gfunci, dreal ww,
                              matfunc_dcmplx *GxGfunc, matfunc_dcmplx *GxGfunc_aug)
{
  int          ngrd, ig, ngrd_aug, ig_aug, ntag, nrow, ncol, isub;
  dreal        dxval, dxsub, xsub, xsubmin;
  mat2d_dcmplx *GxGmat = NULL, *GxGmat_aug = NULL;
  
  ngrd = GxGfunc->ngrd;
  ntag = GxGfunc->ntag;
  ngrd_aug = ngrd + (NSUB - 1) * ntag;
  
  nrow = GxGfunc->nrow;
  ncol = GxGfunc->ncol;
  
  GxGfunc_aug->deall(GxGfunc_aug);
  GxGfunc_aug->alloc(ngrd_aug, nrow, ncol, GxGfunc_aug);
  
  for(ig = 0, ig_aug = 0; ig < ngrd - 1; ++ig, ++ig_aug) {
    GxGfunc_aug->xval[ig_aug] = GxGfunc->xval[ig];
    if(GxGfunc->tag[ig]) {
      dxval = GxGfunc->xval[ig+1] - GxGfunc->xval[ig];
      dxsub = dxval / NSUB;
      xsubmin = GxGfunc->xval[ig];
      for(isub = 1; isub < NSUB; ++isub) {
        ++ig_aug;
        xsub = xsubmin + dxsub * isub;
        GxGfunc_aug->xval[ig_aug] = xsub;
      }
    }
  }
  
  for(ig = 0, ig_aug = 0; ig < ngrd - 1; ++ig, ++ig_aug) {
    GxGmat = GxGfunc->ptr[ig];
    GxGmat_aug = GxGfunc_aug->ptr[ig_aug];
    
    mat2d_copy_dcmplx(GxGmat, GxGmat_aug);
    if(GxGfunc->tag[ig]) {
      for(isub = 1; isub < NSUB; ++isub) {
        ++ig_aug;
        xsub = GxGfunc_aug->xval[ig_aug];
        GxGmat_aug = GxGfunc_aug->ptr[ig_aug];
        nca_gxgmat(Gfuncj, Gfunci, ww, xsub, GxGmat_aug);
      }
    }
  }
  
  GxGfunc_aug->xval[ngrd_aug-1] = GxGfunc->xval[ngrd-1];
  mat2d_copy_dcmplx(GxGfunc->ptr[ngrd-1], GxGfunc_aug->ptr[ngrd_aug-1]);
}

static matfunc_dcmplx *nca_gxgfunc_adpt(int jpspblkj, int ipspblkj, matfunc_dcmplx *Grnpspblkj,
                                        int jpspblki, int ipspblki, matfunc_dcmplx *Grnpspblki,
                                        dreal ww)
{
  int            ngrdj, ngrdi, ig, ivib, jvib, cnvg;
  mat2d_dcmplx   *Grnpspmat = NULL, *Gmat = NULL;
  matfunc_dcmplx *Gfuncj = NULL, *Gfunci = NULL,
                 *GxGfunc = NULL, *GxGfunc_aug = NULL;
  
  // Gfuncj, subblock index jpspblkj, ipspblkj
  ngrdj = Grnpspblkj->ngrd;
  Gfuncj = matfunc_new_dcmplx();
  check_mem(Gfuncj, "Gfuncj");
  Gfuncj->alloc(ngrdj, NVIB, NVIB, Gfuncj);
  for(ig = 0; ig < ngrdj; ++ig) {
    Gfuncj->xval[ig] = Grnpspblkj->xval[ig];
    Gmat = Gfuncj->ptr[ig];
    Grnpspmat = Grnpspblkj->ptr[ig];
    for(ivib = 0; ivib < NVIB; ++ivib) {
      for(jvib = 0; jvib < NVIB; ++jvib) {
        Gmat->ptr[ivib][jvib] = Grnpspmat->ptr[ipspblkj*NVIB+ivib][jpspblkj*NVIB+jvib];
      }
    }
  }
  
  // Gfunci, subblock index jpspblki, ipspblki
  ngrdi = Grnpspblki->ngrd;
  Gfunci = matfunc_new_dcmplx();
  check_mem(Gfunci, "Gfunci");
  Gfunci->alloc(ngrdi, NVIB, NVIB, Gfunci);
  for(ig = 0; ig < ngrdi; ++ig) {
    Gfunci->xval[ig] = Grnpspblki->xval[ig];
    Gmat = Gfunci->ptr[ig];
    Grnpspmat = Grnpspblki->ptr[ig];
    for(ivib = 0; ivib < NVIB; ++ivib) {
      for(jvib = 0; jvib < NVIB; ++jvib) {
        Gmat->ptr[ivib][jvib] = Grnpspmat->ptr[ipspblki*NVIB+ivib][jpspblki*NVIB+jvib];
      }
    }
  }
  
  GxGfunc = matfunc_new_dcmplx();
  
  nca_gxgfunc_unifm(Gfuncj, Gfunci, ww, GxGfunc);
  
  for(cnvg = 0; GxGfunc->ngrd < MAXNGRD && !cnvg; ) {
    //nca_intrp_tag(GxGfunc);
    matfunc_intrp_tag_dcmplx(GxGfunc);
    if(GxGfunc->ntag) {
      GxGfunc_aug = matfunc_new_dcmplx(); check_mem(GxGfunc_aug, "GxGfunc_aug");
      nca_gxgfunc_intrp(Gfuncj, Gfunci, ww, GxGfunc, GxGfunc_aug);
      matfunc_del_dcmplx(GxGfunc);
      GxGfunc = GxGfunc_aug;
    } else {
      cnvg = 1;
    }
  }
  
  matfunc_del_dcmplx(Gfunci);
  matfunc_del_dcmplx(Gfuncj);
  
  return GxGfunc;
  
 error:
  return NULL;
}

static void nca_grnlqusmat(int ispin, dreal ww, mat2d_dcmplx *GrnLqusmat)
{
  const int      thr = 1e-9;
  int            iq, jq, iblk, jblk, npspblki, npspblkj, ipspblki, jpspblki, ipspblkj, jpspblkj,
                 ipspi, jpspi, ipspj, jpspj, ivib, jvib, ngrd, ig;
  dcmplx         xil, xir;
  dcmplx         *gxg = NULL;
  matfunc_dcmplx *GrnLpspblk = NULL, *GrnGpspblk = NULL;
  matfunc_dcmplx *GxGfunc = NULL;
  
  GrnLqusmat->reset(GrnLqusmat);
  
  for(iq = 0; iq < Nqus; ++iq) {
    for(jq = 0; jq < Nqus; ++jq) {
      
      for(iblk = 0; iblk < Nblk; ++iblk) {
        npspblki   = Npspblk[iblk];
        GrnGpspblk = GrnGpsp_aug[iblk]; // use augmented version
        for(ipspblki = 0; ipspblki < npspblki; ++ipspblki) {
          ipspi = Pspblk[iblk][ipspblki];
          for(jpspblki = 0; jpspblki < npspblki; ++jpspblki) {
            jpspi = Pspblk[iblk][jpspblki];
            
            for(jblk = 0; jblk < Nblk; ++jblk) {
              npspblkj   = Npspblk[jblk];
              GrnLpspblk = GrnLpsp_aug[jblk]; // use augmented version
              for(ipspblkj = 0; ipspblkj < npspblkj; ++ipspblkj) {
                ipspj = Pspblk[jblk][ipspblkj];
                for(jpspblkj = 0; jpspblkj < npspblkj; ++jpspblkj) {
                  jpspj = Pspblk[jblk][jpspblkj];
                  
                  xil = Xipsp[ispin][iq]->ptr[jpspi][ipspj];
                  xir = Xipsp[ispin][jq]->ptr[ipspi][jpspj];
                  
                  if(cabs(xil * conj(xir)) > thr) {
                    GxGfunc = nca_gxgfunc_adpt(jpspblkj, ipspblkj, GrnLpspblk,
                                               jpspblki, ipspblki, GrnGpspblk, ww);
                    ngrd = GxGfunc->ngrd;
                    gxg = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
                    for(ivib = 0; ivib < NVIB; ++ivib) {
                      for(jvib = 0; jvib < NVIB; ++jvib) {
                        for(ig = 0; ig < ngrd; ++ig) {
                          gxg[ig] = GxGfunc->ptr[ig]->ptr[ivib][jvib];
                        }
                        GrnLqusmat->ptr[iq][jq] +=-I / (2.0*PI) * Zeta[jpspj] * xil * conj(xir)
                          * numeric_quadrat_dcmplx(ngrd, GxGfunc->xval, gxg);
                      }
                    }
                    freeup(gxg);
                    matfunc_del_dcmplx(GxGfunc);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

static void nca_grngqusmat(int ispin, dreal ww, mat2d_dcmplx *GrnGqusmat)
{
  const int      thr = 1e-9;
  int            iq, jq, iblk, jblk, npspblki, npspblkj, ipspblki, jpspblki, ipspblkj, jpspblkj,
                 ipspi, jpspi, ipspj, jpspj, ivib, jvib, ngrd, ig;
  dcmplx         xil, xir;
  dcmplx         *gxg = NULL;
  matfunc_dcmplx *GrnGpspblk = NULL, *GrnLpspblk = NULL;
  matfunc_dcmplx *GxGfunc = NULL;
  
  GrnGqusmat->reset(GrnGqusmat);
  
  for(iq = 0; iq < Nqus; ++iq) {
    for(jq = 0; jq < Nqus; ++jq) {
      
      for(iblk = 0; iblk < Nblk; ++iblk) {
        npspblki   = Npspblk[iblk];
        GrnLpspblk = GrnLpsp_aug[iblk]; // use augmented version
        for(ipspblki = 0; ipspblki < npspblki; ++ipspblki) {
          ipspi = Pspblk[iblk][ipspblki];
          for(jpspblki = 0; jpspblki < npspblki; ++jpspblki) {
            jpspi = Pspblk[iblk][jpspblki];
            
            for(jblk = 0; jblk < Nblk; ++jblk) {
              npspblkj   = Npspblk[jblk];
              GrnGpspblk = GrnGpsp_aug[jblk]; // use augmented version
              for(ipspblkj = 0; ipspblkj < npspblkj; ++ipspblkj) {
                ipspj = Pspblk[jblk][ipspblkj];
                for(jpspblkj = 0; jpspblkj < npspblkj; ++jpspblkj) {
                  jpspj = Pspblk[jblk][jpspblkj];
                  
                  xil = Xipsp[ispin][iq]->ptr[jpspi][ipspj];
                  xir = Xipsp[ispin][jq]->ptr[ipspi][jpspj];
                  
                  if(cabs(xil * conj(xir)) > thr) {
                    GxGfunc = nca_gxgfunc_adpt(jpspblkj, ipspblkj, GrnGpspblk,
                                               jpspblki, ipspblki, GrnLpspblk, ww);
                    ngrd = GxGfunc->ngrd;
                    gxg = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
                    for(ivib = 0; ivib < NVIB; ++ivib) {
                      for(jvib = 0; jvib < NVIB; ++jvib) {
                        for(ig = 0; ig < ngrd; ++ig) {
                          gxg[ig] = GxGfunc->ptr[ig]->ptr[ivib][jvib];
                        }
                        GrnGqusmat->ptr[iq][jq] +=-I / (2.0*PI) * Zeta[jpspj] * xil * conj(xir)
                          * numeric_quadrat_dcmplx(ngrd, GxGfunc->xval, gxg);
                      }
                    }
                    freeup(gxg);
                    matfunc_del_dcmplx(GxGfunc);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

static matfunc_dcmplx *nca_grnqus_unifm(int ispin, gqsmatfunc nca_grnqusmat)
{
  int            ig;
  size_t         nlen;
  dreal          dE, Eng;
  dcmplx         *Grnbuf = NULL;
  mat2d_dcmplx   *Grnqusmat = NULL;
  matfunc_dcmplx *Grnqusfunc = NULL;
  
  nlen = sizeof(dcmplx)*Nqus*Nqus;
  
  Grnbuf = (dcmplx *) calloc(Nqus*Nqus*Nmesh, sizeof(dcmplx));
  
  Grnqusfunc = matfunc_new_dcmplx();
  check_mem(Grnqusfunc, "Grnqusfunc");
  Grnqusfunc->alloc(Nmesh, Nqus, Nqus, Grnqusfunc);
  
  dE = (Emax - Emin) / (Nmesh - 1);
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    Grnqusfunc->xval[ig] = Eng;
  }
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    Grnqusmat = Grnqusfunc->ptr[ig];
    if(ig % Size == Rank) {
      nca_grnqusmat(ispin, Eng, Grnqusmat);
      memcpy(Grnbuf + Nqus*Nqus*ig, Grnqusmat->addr, nlen);
    }
  }
  
  sys_allreduce_sum_dcmplx(Grnbuf, Nqus*Nqus*Nmesh); sys_sync();
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Grnqusmat = Grnqusfunc->ptr[ig];
    if(ig % Size != Rank) {
      memcpy(Grnqusmat->addr, Grnbuf + Nqus*Nqus*ig, nlen);
    }
  }
  
  freeup(Grnbuf);
  
  return Grnqusfunc;
  
 error:
  return NULL;
}

static matfunc_dcmplx *nca_grnqus_ret(matfunc_dcmplx *GrnGqusfunc,
                                      matfunc_dcmplx *GrnLqusfunc)
{
  int            ig, iq, jq;
  dreal          dE, Eng;
  dcmplx         *dgrn = NULL, *grnr = NULL;
  mat2d_dcmplx   *GrnGqusmat = NULL, *GrnLqusmat = NULL,
                 *GrnRqusmat = NULL; 
  matfunc_dcmplx *GrnRqusfunc = NULL;
  
  GrnRqusfunc = matfunc_new_dcmplx();
  check_mem(GrnRqusfunc, "GrnRqusfunc");
  GrnRqusfunc->alloc(Nmesh, Nqus, Nqus, GrnRqusfunc);
  
  dE = (Emax - Emin) / (Nmesh - 1);
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    GrnRqusfunc->xval[ig] = Eng;
  }
  
  dgrn = (dcmplx *) calloc(Nmesh_aug-1, sizeof(dcmplx));
  grnr = (dcmplx *) calloc(Nmesh_aug-1, sizeof(dcmplx));
  
  for(iq = 0 ; iq < Nqus; ++iq) {
    for(jq = 0; jq < Nqus; ++jq) {
      for(ig = 0; ig < Nmesh; ++ig) {
        GrnGqusmat = GrnGqusfunc->ptr[ig];
        GrnLqusmat = GrnLqusfunc->ptr[ig];
        dgrn[ig+(Nmesh-1)/2] = GrnGqusmat->ptr[iq][jq] - GrnLqusmat->ptr[iq][jq];
      }
      numeric_kronig(Nmesh_aug-1, dgrn, grnr);
      for(ig = 0; ig < Nmesh; ++ig) {
        GrnRqusmat = GrnRqusfunc->ptr[ig];
        GrnRqusmat->ptr[iq][jq] = grnr[ig+(Nmesh-1)/2];
      }
    }
  }
  
  freeup(grnr); freeup(dgrn);
  
  return GrnRqusfunc;
  
 error:
  return NULL;
}

static matfunc_dcmplx *nca_grnqus_adv(matfunc_dcmplx *GrnRqusfunc)
{
  int            ig, iq, jq;
  dreal          dE, Eng;
  mat2d_dcmplx   *GrnRqusmat = NULL, *GrnAqusmat = NULL; 
  matfunc_dcmplx *GrnAqusfunc = NULL;
  
  GrnAqusfunc = matfunc_new_dcmplx();
  check_mem(GrnAqusfunc, "GrnAqusfunc");
  GrnAqusfunc->alloc(Nmesh, Nqus, Nqus, GrnAqusfunc);
  
  dE = (Emax - Emin) / (Nmesh - 1);
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    GrnAqusfunc->xval[ig] = Eng;
  }
  
  for(ig = 0; ig < Nmesh; ++ig) {
    for(iq = 0 ; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        GrnRqusmat = GrnRqusfunc->ptr[ig];
        GrnAqusmat = GrnAqusfunc->ptr[ig];
        GrnAqusmat->ptr[iq][jq] = conj(GrnRqusmat->ptr[jq][iq]);
      }
    }
  }
  
  return GrnAqusfunc;
  
 error:
  return NULL;
}

static void nca_dosqus_out(void)
{
  FILE  *fp = NULL;
  int   ispin, iq, ngrd, ig;
  char  *filename = NULL;
  dreal *dos = NULL;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      ngrd = GrnRqus[ispin]->ngrd;
      dos = (dreal *) calloc(ngrd, sizeof(dreal));
      for(ig = 0; ig < ngrd; ++ig) {
        dos[ig] =-cimag(GrnRqus[ispin]->ptr[ig]->ptr[iq][iq]) / PI;
      }
      // Output
      filename = (char *) calloc(SHRT_MAX, sizeof(char));
      sprintf(filename, "%s_%1d_%04d.ascii", "DosqusNCA", ispin, iq);
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", GrnRqus[ispin]->xval[ig], dos[ig]);
        }
        fclose(fp);
      }
      freeup(filename);
      freeup(dos);
    }  
  }
}

static void nca_popqus_out(void)
{
  FILE  *fp = NULL;
  int   ispin, iq, ngrd, ig;
  char  *filename = NULL;
  dreal *pop = NULL;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      ngrd = GrnLqus[ispin]->ngrd;
      pop = (dreal *) calloc(ngrd, sizeof(dreal));
      for(ig = 0; ig < ngrd; ++ig) {
        pop[ig] = cimag(GrnLqus[ispin]->ptr[ig]->ptr[iq][iq]) / (2.0*PI);
      }
      // Output
      filename = (char *) calloc(SHRT_MAX, sizeof(char));
      sprintf(filename, "%s_%1d_%04d.ascii", "PopqusNCA", ispin, iq);
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", GrnLqus[ispin]->xval[ig], pop[ig]);
        }
        fclose(fp);
      }
      freeup(filename);
      freeup(pop);
    }  
  }
}

void nca_grnqus(void)
{
  int ispin;
  
  nca_psp_aug_alloc();
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    matfunc_del_dcmplx(GrnLqus[ispin]); GrnLqus[ispin] = nca_grnqus_unifm(ispin, nca_grnlqusmat);
    matfunc_del_dcmplx(GrnGqus[ispin]); GrnGqus[ispin] = nca_grnqus_unifm(ispin, nca_grngqusmat);
    matfunc_del_dcmplx(GrnRqus[ispin]); GrnRqus[ispin] = nca_grnqus_ret(GrnGqus[ispin], GrnLqus[ispin]);
    matfunc_del_dcmplx(GrnAqus[ispin]); GrnAqus[ispin] = nca_grnqus_adv(GrnRqus[ispin]);
  }
  
  nca_dosqus_out();
  nca_popqus_out();
  
  nca_psp_aug_deall();
}
