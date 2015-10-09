#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "dbg.h"
#include "sys.h"
#include "lapack.h"
#include "numeric.h"
#include "constants.h"
#include "variables.h"
#include "stats.h"
#include "perturb.h"
#include "nca.h"

matfunc_dcmplx **PiR_aug = NULL, **PiA_aug = NULL, **PiL_aug = NULL, **PiG_aug = NULL,
               **LamR_aug = NULL, **LamA_aug = NULL, **LamL_aug = NULL, **LamG_aug = NULL;

static const int      MAXITER = 20;
static const dreal    THRSHLD = 1e-6;
static const dreal    TOLDPOP = 2e-9;
static dreal          **SumPopqus_old = NULL, **SumPopqus_new = NULL;
static matfunc_dcmplx **GrnM = NULL;

static const char *dmftlabel = "dmft:           ";

void dmft_matfunc_aug(const matfunc_dcmplx *matfunc, matfunc_dcmplx *matfunc_aug)
{
  int    nrow, ncol, ngrd, ngrd_aug, ngrd_neg, ngrd_pos, ig_aug, ig;
  size_t nlen;
  dreal  *xval = NULL;
  
  ngrd = matfunc->ngrd;
  ngrd_aug = NAUG * (ngrd - 1) + 1;
  
  nrow = matfunc->nrow;
  ncol = matfunc->ncol;
  
  nlen = sizeof(dcmplx)*nrow*ncol;
  
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
    memset(matfunc_aug->ptr[ig_aug]->addr, 0, nlen);
  }
  
  for(ig_aug = ngrd_pos; ig_aug < ngrd_pos+ngrd; ++ig_aug) {
    ig = ig_aug - ngrd_pos;
    matfunc_aug->xval[ig_aug] = matfunc->xval[ig];
    mat2d_copy_dcmplx(matfunc->ptr[ig], matfunc_aug->ptr[ig_aug]);
  }
  
  for(ig_aug = ngrd_pos+ngrd; ig_aug < ngrd_aug; ++ig_aug) {
    ig = ig_aug - ngrd_pos - ngrd + 1;
    matfunc_aug->xval[ig_aug] = matfunc->xval[ig] + (Emax - Emin);
    memset(matfunc_aug->ptr[ig_aug]->addr, 0, nlen);
  }
  
  return;
  
 error:
  abort();
}

void dmft_alloc(void)
{
  int ispin;
  // Population, old
  SumPopqus_old = (dreal **) calloc(Nspin, sizeof(dreal *));
  SumPopqus_new = (dreal **) calloc(Nspin, sizeof(dreal *));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    SumPopqus_old[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
    SumPopqus_new[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
  }
  // Matrix form of Green's functions
  GrnM = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnM[ispin] = matfunc_new_dcmplx();
    check_mem(GrnM[ispin], "GrnM[ispin]");
    GrnM[ispin]->alloc(Nmesh, NDIMSUB*Ndim, NDIMSUB*Ndim, GrnM[ispin]);
  }
  // Augmented hybridization
  LamR_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  LamA_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  LamL_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  LamG_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    LamR_aug[ispin] = matfunc_new_dcmplx();
    check_mem(LamR_aug[ispin], "LamR_aug[ispin]");
    LamR_aug[ispin]->alloc(Nmesh_aug, Nqus, Nqus, LamR_aug[ispin]);
    
    LamA_aug[ispin] = matfunc_new_dcmplx();
    check_mem(LamA_aug[ispin], "LamA_aug[ispin]");
    LamA_aug[ispin]->alloc(Nmesh_aug, Nqus, Nqus, LamA_aug[ispin]);
    
    LamL_aug[ispin] = matfunc_new_dcmplx();
    check_mem(LamL_aug[ispin], "LamL_aug[ispin]");
    LamL_aug[ispin]->alloc(Nmesh_aug, Nqus, Nqus, LamL_aug[ispin]);
    
    LamG_aug[ispin] = matfunc_new_dcmplx();
    check_mem(LamG_aug[ispin], "LamG_aug[ispin]");
    LamG_aug[ispin]->alloc(Nmesh_aug, Nqus, Nqus, LamG_aug[ispin]);
  }
  // Augmented optic coupling
  if(isOptic) {
    PiR_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    PiA_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    PiL_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    PiG_aug = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    for(ispin = 0; ispin < Nspin; ++ispin) {
      PiR_aug[ispin] = matfunc_new_dcmplx();
      check_mem(PiR_aug[ispin], "PiR_aug[ispin]");
      PiR_aug[ispin]->alloc(Nmesh_aug, 1, 1, PiR_aug[ispin]);
      
      PiA_aug[ispin] = matfunc_new_dcmplx();
      check_mem(PiA_aug[ispin], "PiA_aug[ispin]");
      PiA_aug[ispin]->alloc(Nmesh_aug, 1, 1, PiA_aug[ispin]);
      
      PiL_aug[ispin] = matfunc_new_dcmplx();
      check_mem(PiL_aug[ispin], "PiL_aug[ispin]");
      PiL_aug[ispin]->alloc(Nmesh_aug, 1, 1, PiL_aug[ispin]);
      
      PiG_aug[ispin] = matfunc_new_dcmplx();
      check_mem(PiG_aug[ispin], "PiG_aug[ispin]");
      PiG_aug[ispin]->alloc(Nmesh_aug, 1, 1, PiG_aug[ispin]);
    }
  }
  // Quasiparticle Green's functions and self energies
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnRqus[ispin] = matfunc_new_dcmplx();
    check_mem(GrnRqus[ispin], "GrnRqus[ispin]");
    GrnRqus[ispin]->alloc(Nmesh, Nqus, Nqus, GrnRqus[ispin]);
    
    GrnAqus[ispin] = matfunc_new_dcmplx();
    check_mem(GrnAqus[ispin], "GrnAqus[ispin]");
    GrnAqus[ispin]->alloc(Nmesh, Nqus, Nqus, GrnAqus[ispin]);
    
    GrnLqus[ispin] = matfunc_new_dcmplx();
    check_mem(GrnLqus[ispin], "GrnLqus[ispin]");
    GrnLqus[ispin]->alloc(Nmesh, Nqus, Nqus, GrnLqus[ispin]);
    
    GrnGqus[ispin] = matfunc_new_dcmplx();
    check_mem(GrnGqus[ispin], "GrnGqus[ispin]");
    GrnGqus[ispin]->alloc(Nmesh, Nqus, Nqus, GrnGqus[ispin]);
  }
  for(ispin = 0; ispin < Nspin; ++ispin) { 
    SgmRqus[ispin] = matfunc_new_dcmplx();
    check_mem(SgmRqus[ispin], "SgmRqus[ispin]");
    SgmRqus[ispin]->alloc(Nmesh, Nqus, Nqus, SgmRqus[ispin]);
    
    SgmAqus[ispin] = matfunc_new_dcmplx();
    check_mem(SgmAqus[ispin], "SgmAqus[ispin]");
    SgmAqus[ispin]->alloc(Nmesh, Nqus, Nqus, SgmAqus[ispin]);
    
    SgmLqus[ispin] = matfunc_new_dcmplx();
    check_mem(SgmLqus[ispin], "SgmLqus[ispin]");
    SgmLqus[ispin]->alloc(Nmesh, Nqus, Nqus, SgmLqus[ispin]);
    
    SgmGqus[ispin] = matfunc_new_dcmplx();
    check_mem(SgmGqus[ispin], "SgmGqus[ispin]");
    SgmGqus[ispin]->alloc(Nmesh, Nqus, Nqus, SgmGqus[ispin]);
  }
  // Hybridization potential for quasiparticles
  for(ispin = 0; ispin < Nspin; ++ispin) {
    LamR[ispin] = matfunc_new_dcmplx();
    check_mem(LamR[ispin], "LamR[ispin]");
    LamR[ispin]->alloc(Nmesh, Nqus, Nqus, LamR[ispin]);
    
    LamA[ispin] = matfunc_new_dcmplx();
    check_mem(LamA[ispin], "LamA[ispin]");
    LamA[ispin]->alloc(Nmesh, Nqus, Nqus, LamA[ispin]);
    
    LamL[ispin] = matfunc_new_dcmplx();
    check_mem(LamL[ispin], "LamL[ispin]");
    LamL[ispin]->alloc(Nmesh, Nqus, Nqus, LamL[ispin]);
    
    LamG[ispin] = matfunc_new_dcmplx();
    check_mem(LamG[ispin], "LamG[ispin]");
    LamG[ispin]->alloc(Nmesh, Nqus, Nqus, LamG[ispin]);
  }
  // Optic coupling for quasiparticles
  if(isOptic) {
    for(ispin = 0; ispin < Nspin; ++ispin) {
      PiR[ispin] = matfunc_new_dcmplx();
      check_mem(PiR[ispin], "PiR[ispin]");
      PiR[ispin]->alloc(Nmesh, 1, 1, PiR[ispin]);
      
      PiA[ispin] = matfunc_new_dcmplx();
      check_mem(PiA[ispin], "PiA[ispin]");
      PiA[ispin]->alloc(Nmesh, 1, 1, PiA[ispin]);
      
      PiL[ispin] = matfunc_new_dcmplx();
      check_mem(PiL[ispin], "PiL[ispin]");
      PiL[ispin]->alloc(Nmesh, 1, 1, PiL[ispin]);
      
      PiG[ispin] = matfunc_new_dcmplx();
      check_mem(PiG[ispin], "PiG[ispin]");
      PiG[ispin]->alloc(Nmesh, 1, 1, PiG[ispin]);
    }
  }
  
  return;
  
 error:
  abort();
}

void dmft_deall(void)
{
  int ispin;
  
  if(isOptic) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(PiG[ispin]);
      matfunc_del_dcmplx(PiL[ispin]);
      matfunc_del_dcmplx(PiA[ispin]);
      matfunc_del_dcmplx(PiR[ispin]);
    }
  }
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    matfunc_del_dcmplx(LamG[ispin]);
    matfunc_del_dcmplx(LamL[ispin]);
    matfunc_del_dcmplx(LamA[ispin]);
    matfunc_del_dcmplx(LamR[ispin]);
  }
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    matfunc_del_dcmplx(SgmGqus[ispin]);
    matfunc_del_dcmplx(SgmLqus[ispin]);
    matfunc_del_dcmplx(SgmAqus[ispin]);
    matfunc_del_dcmplx(SgmRqus[ispin]);
  }
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    matfunc_del_dcmplx(GrnGqus[ispin]);
    matfunc_del_dcmplx(GrnLqus[ispin]);
    matfunc_del_dcmplx(GrnAqus[ispin]);
    matfunc_del_dcmplx(GrnRqus[ispin]);
  }
  
  if(isOptic) {
    if(PiG_aug) {
      for(ispin = Nspin - 1; ispin > -1; --ispin) {
        matfunc_del_dcmplx(PiG_aug[ispin]);
      }
      freeup(PiG_aug);
    }
    if(PiL_aug) {
      for(ispin = Nspin - 1; ispin > -1; --ispin) {
        matfunc_del_dcmplx(PiL_aug[ispin]);
      }
      freeup(PiL_aug);
    }
    if(PiA_aug) {
      for(ispin = Nspin - 1; ispin > -1; --ispin) {
        matfunc_del_dcmplx(PiA_aug[ispin]);
      }
      freeup(PiA_aug);
    }
    if(PiR_aug) {
      for(ispin = Nspin - 1; ispin > -1; --ispin) {
        matfunc_del_dcmplx(PiR_aug[ispin]);
      }
      freeup(PiR_aug);
    }
  }
  
  if(LamG_aug) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamG_aug[ispin]);
    }
    freeup(LamG_aug);
  }
  if(LamL_aug) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamL_aug[ispin]);
    }
    freeup(LamL_aug);
  }
  if(LamA_aug) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamA_aug[ispin]);
    }
    freeup(LamA_aug);
  }
  if(LamR_aug) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamR_aug[ispin]);
    }
    freeup(LamR_aug);
  }
  
  if(GrnM) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnM[ispin]);
    }
    freeup(GrnM);
  }
  
  if(SumPopqus_new) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      if(SumPopqus_new[ispin]) freeup(SumPopqus_new[ispin]);
    }
    freeup(SumPopqus_new);
  }
  if(SumPopqus_old) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      if(SumPopqus_old[ispin]) freeup(SumPopqus_old[ispin]);
    }
    freeup(SumPopqus_old);
  }
}

static void dmft_sgmrmat(const mat2d_dcmplx *SgmRqusmat, mat2d_dcmplx *SgmRmat)
{
  int    idim, jdim, ienv, iq, jq;
  dreal  gm;
  dcmplx se;
  
  SgmRmat->reset(SgmRmat);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(ienv = 0; ienv < Nenv; ++ienv) {
      gm = Gamma[ienv][idim];
      se =-0.5 * I * gm;
      SgmRmat->ptr[idim][idim] += se;
    }
    SgmRmat->ptr[idim][idim] +=-I * INFTES;
  }
  
  for(iq = 0; iq < Nqus; ++iq) {
    idim = Indqus[iq];
    for(jq = 0; jq < Nqus; ++jq) {
      jdim = Indqus[jq];
      SgmRmat->ptr[idim][jdim] += SgmRqusmat->ptr[iq][jq];
    }
  }
}

static void dmft_sgmamat(const mat2d_dcmplx *SgmAqusmat, mat2d_dcmplx *SgmAmat)
{
  int    idim, jdim, ienv, iq, jq;
  dreal  gm;
  dcmplx se;
  
  SgmAmat->reset(SgmAmat);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(ienv = 0; ienv < Nenv; ++ienv) {
      gm = Gamma[ienv][idim];
      se = 0.5 * I * gm;
      SgmAmat->ptr[idim][idim] += se;
    }
    SgmAmat->ptr[idim][idim] += I * INFTES;
  }
  
  for(iq = 0; iq < Nqus; ++iq) {
    idim = Indqus[iq];
    for(jq = 0; jq < Nqus; ++jq) {
      jdim = Indqus[jq];
      SgmAmat->ptr[idim][jdim] += SgmAqusmat->ptr[iq][jq];
    }
  }
}

static void dmft_sgmlmat(dreal ww, const mat2d_dcmplx *SgmLqusmat, mat2d_dcmplx *SgmLmat)
{
  int    idim, jdim, ienv, iq, jq;
  dreal  gm, mu, ff;
  dcmplx se;
  
  SgmLmat->reset(SgmLmat);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(ienv = 0; ienv < Nenv; ++ienv) {
      gm = Gamma[ienv][idim];
      mu = Muenv[ienv];
      ff = stats_ferm(Beta, mu, ww);
      se = I * ff * gm;
      SgmLmat->ptr[idim][idim] += se;
    }
    SgmLmat->ptr[idim][idim] += I * INFTES;
  }
  
  for(iq = 0; iq < Nqus; ++iq) {
    idim = Indqus[iq];
    for(jq = 0; jq < Nqus; ++jq) {
      jdim = Indqus[jq];
      SgmLmat->ptr[idim][jdim] += SgmLqusmat->ptr[iq][jq];
    }
  }
}

static void dmft_sgmgmat(dreal ww, const mat2d_dcmplx *SgmGqusmat, mat2d_dcmplx *SgmGmat)
{
  int    idim, jdim, ienv, iq, jq;
  dreal  gm, mu, ff;
  dcmplx se;
  
  SgmGmat->reset(SgmGmat);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(ienv = 0; ienv < Nenv; ++ienv) {
      gm = Gamma[ienv][idim];
      mu = Muenv[ienv];
      ff = stats_invsferm(Beta, mu, ww);
      se =-I * ff * gm;
      SgmGmat->ptr[idim][idim] += se;
    }
    SgmGmat->ptr[idim][idim] +=-I * INFTES;
  }
  
  for(iq = 0; iq < Nqus; ++iq) {
    idim = Indqus[iq];
    for(jq = 0; jq < Nqus; ++jq) {
      jdim = Indqus[jq];
      SgmGmat->ptr[idim][jdim] += SgmGqusmat->ptr[iq][jq];
    }
  }
}

static void dmft_grnm(int ispin, matfunc_dcmplx *GrnMfunc)
{
  int          ig, idim, jdim, idimsub, jdimsub;
  dreal        dE, Eng;
  mat2d_dcmplx *SgmRqusmat = NULL, *SgmAqusmat = NULL,
               *SgmLqusmat = NULL, *SgmGqusmat = NULL;
  mat2d_dcmplx *SgmRmat = NULL, *SgmAmat = NULL,
               *SgmLmat = NULL, *SgmGmat = NULL,
               *GrnMmat = NULL;
  
  SgmRmat = mat2d_new_dcmplx();
  check_mem(SgmRmat, "SgmRmat");
  SgmRmat->alloc(Ndim, Ndim, SgmRmat);
  
  SgmAmat = mat2d_new_dcmplx();
  check_mem(SgmAmat, "SgmAmat");
  SgmAmat->alloc(Ndim, Ndim, SgmAmat);
  
  SgmLmat = mat2d_new_dcmplx();
  check_mem(SgmLmat, "SgmLmat");
  SgmLmat->alloc(Ndim, Ndim, SgmLmat);
  
  SgmGmat = mat2d_new_dcmplx();
  check_mem(SgmGmat, "SgmGmat");
  SgmGmat->alloc(Ndim, Ndim, SgmGmat);
  
  dE = (Emax - Emin) / (Nmesh - 1);
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin+ dE * ig;
    GrnMfunc->xval[ig] = Eng;
    GrnMmat = GrnMfunc->ptr[ig];
    
    SgmRqusmat = SgmRqus[ispin]->ptr[ig];
    SgmAqusmat = SgmAqus[ispin]->ptr[ig];
    SgmLqusmat = SgmLqus[ispin]->ptr[ig];
    SgmGqusmat = SgmGqus[ispin]->ptr[ig];
    
    dmft_sgmrmat(SgmRqusmat, SgmRmat);
    dmft_sgmamat(SgmAqusmat, SgmAmat);
    dmft_sgmlmat(Eng, SgmLqusmat, SgmLmat);
    dmft_sgmgmat(Eng, SgmGqusmat, SgmGmat);
    // Loop over subblock
    for(idim = 0; idim < Ndim; ++idim) {
      for(jdim = 0; jdim < Ndim; ++jdim) {
        // Part 1
        jdimsub = NDIMSUB * jdim;
        idimsub = NDIMSUB * idim;
        GrnMmat->ptr[idimsub][jdimsub] =-Hmat->ptr[idim][jdim]
          - SgmLmat->ptr[idim][jdim] - SgmRmat->ptr[idim][jdim];
        // Part 2
        jdimsub = NDIMSUB * jdim;
        idimsub = NDIMSUB * idim + 1;
        GrnMmat->ptr[idimsub][jdimsub] = SgmLmat->ptr[idim][jdim];
        // Part 3
        jdimsub = NDIMSUB * jdim + 1;
        idimsub = NDIMSUB * idim;
        GrnMmat->ptr[idimsub][jdimsub] = SgmGmat->ptr[idim][jdim];
        // Part 4
        jdimsub = NDIMSUB * jdim + 1;
        idimsub = NDIMSUB * idim + 1;
        GrnMmat->ptr[idimsub][jdimsub] = Hmat->ptr[idim][jdim]
          - SgmLmat->ptr[idim][jdim] + SgmAmat->ptr[idim][jdim];
      }
      // Diagonal part
      idimsub = NDIMSUB * idim;
      GrnMmat->ptr[idimsub][idimsub] += Eng;
      idimsub = NDIMSUB * idim + 1;
      GrnMmat->ptr[idimsub][idimsub] -= Eng;
    }
    lapack_zinvs(NDIMSUB*Ndim, GrnMmat->addr);
  }
  
  mat2d_del_dcmplx(SgmGmat); mat2d_del_dcmplx(SgmLmat);
  mat2d_del_dcmplx(SgmAmat); mat2d_del_dcmplx(SgmRmat);
  
  return;
  
 error:
  abort();
}

static void dmft_grnqus(const matfunc_dcmplx *GrnMfunc, 
                        matfunc_dcmplx *GrnRqusfunc, matfunc_dcmplx *GrnAqusfunc,
                        matfunc_dcmplx *GrnLqusfunc, matfunc_dcmplx *GrnGqusfunc)
{
  int          ig, iq, jq, idim, jdim, idimsub, jdimsub;
  dreal        *xval = NULL;
  mat2d_dcmplx *GrnRqusmat = NULL, *GrnAqusmat = NULL,
               *GrnLqusmat = NULL, *GrnGqusmat = NULL,
               *GrnMmat = NULL;
  
  xval = GrnMfunc->xval;
  
  for(ig = 0; ig < Nmesh; ++ig) {
    GrnRqusfunc->xval[ig] = xval[ig];
    GrnAqusfunc->xval[ig] = xval[ig];
    GrnLqusfunc->xval[ig] = xval[ig];
    GrnGqusfunc->xval[ig] = xval[ig];
    
    GrnMmat = GrnMfunc->ptr[ig];
    GrnRqusmat = GrnRqusfunc->ptr[ig];
    GrnAqusmat = GrnAqusfunc->ptr[ig];
    GrnLqusmat = GrnLqusfunc->ptr[ig];
    GrnGqusmat = GrnGqusfunc->ptr[ig];
    for(iq = 0; iq < Nqus; ++iq) {
      idim = Indqus[iq];
      idimsub = NDIMSUB * idim;
      for(jq = 0; jq < Nqus; ++jq) {
        jdim = Indqus[jq];
        jdimsub = NDIMSUB * jdim;
        GrnLqusmat->ptr[iq][jq] = GrnMmat->ptr[idimsub+1][jdimsub];
        GrnGqusmat->ptr[iq][jq] = GrnMmat->ptr[idimsub][jdimsub+1];
        GrnRqusmat->ptr[iq][jq] = GrnMmat->ptr[idimsub][jdimsub] - GrnMmat->ptr[idimsub+1][jdimsub];
        GrnAqusmat->ptr[iq][jq] = GrnMmat->ptr[idimsub+1][jdimsub] - GrnMmat->ptr[idimsub+1][jdimsub+1];
      }
    }
  }
}

static void dmft_lam(const matfunc_dcmplx *GrnRqusfunc, const matfunc_dcmplx *GrnAqusfunc,
                     const matfunc_dcmplx *GrnLqusfunc, const matfunc_dcmplx *GrnGqusfunc,
                     const matfunc_dcmplx *SgmRqusfunc, const matfunc_dcmplx *SgmAqusfunc,
                     const matfunc_dcmplx *SgmLqusfunc, const matfunc_dcmplx *SgmGqusfunc,
                     matfunc_dcmplx *LamRfunc, matfunc_dcmplx *LamAfunc,
                     matfunc_dcmplx *LamLfunc, matfunc_dcmplx *LamGfunc)
{
  int          iq, jq, idim, jdim, ig;
  dreal        dE, Eng;
  mat2d_dcmplx *Hqus = NULL, *Grninvs = NULL;
  mat2d_dcmplx *GrnRqusmat = NULL, *GrnAqusmat = NULL,
               *GrnLqusmat = NULL, *GrnGqusmat = NULL,
               *SgmRqusmat = NULL, *SgmAqusmat = NULL,
               *SgmLqusmat = NULL, *SgmGqusmat = NULL,
               *LamRmat = NULL, *LamAmat = NULL,
               *LamLmat = NULL, *LamGmat = NULL;
  
  Hqus = mat2d_new_dcmplx();
  check_mem(Hqus, "Hqus");
  Hqus->alloc(Nqus, Nqus, Hqus);
  
  Grninvs = mat2d_new_dcmplx();
  check_mem(Grninvs, "Grninvs");
  Grninvs->alloc(NDIMSUB*Nqus, NDIMSUB*Nqus, Grninvs);
  
  for(iq = 0; iq < Nqus; ++iq) {
    idim = Indqus[iq];
    for(jq = 0; jq < Nqus; ++jq) {
      jdim = Indqus[jq];
      Hqus->ptr[iq][jq] = Hmat->ptr[idim][jdim];
    }
  }
  
  dE = (Emax - Emin) / (Nmesh - 1);
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    LamRfunc->xval[ig] = Eng; LamAfunc->xval[ig] = Eng;
    LamLfunc->xval[ig] = Eng; LamGfunc->xval[ig] = Eng;
    GrnRqusmat = GrnRqusfunc->ptr[ig]; GrnAqusmat = GrnAqusfunc->ptr[ig];
    GrnLqusmat = GrnLqusfunc->ptr[ig]; GrnGqusmat = GrnGqusfunc->ptr[ig];
    SgmRqusmat = SgmRqusfunc->ptr[ig]; SgmAqusmat = SgmAqusfunc->ptr[ig];
    SgmLqusmat = SgmLqusfunc->ptr[ig]; SgmGqusmat = SgmGqusfunc->ptr[ig];
    LamRmat = LamRfunc->ptr[ig]; LamAmat = LamAfunc->ptr[ig];
    LamLmat = LamLfunc->ptr[ig]; LamGmat = LamGfunc->ptr[ig];
    
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        Grninvs->ptr[iq][jq] = GrnLqusmat->ptr[iq][jq] + GrnRqusmat->ptr[iq][jq];
        Grninvs->ptr[Nqus+iq][jq] = GrnLqusmat->ptr[iq][jq];
        Grninvs->ptr[iq][Nqus+jq] = GrnGqusmat->ptr[iq][jq];
        Grninvs->ptr[Nqus+iq][Nqus+jq] = GrnLqusmat->ptr[iq][jq] - GrnAqusmat->ptr[iq][jq];
      }
    }
    
    lapack_zinvs(NDIMSUB*Nqus, Grninvs->addr);
    
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        LamLmat->ptr[iq][jq] = Grninvs->ptr[Nqus+iq][jq] - SgmLqusmat->ptr[iq][jq];
        LamGmat->ptr[iq][jq] = Grninvs->ptr[iq][Nqus+jq] - SgmGqusmat->ptr[iq][jq];
        LamRmat->ptr[iq][jq] =-Hqus->ptr[iq][jq] - SgmLqusmat->ptr[iq][jq] - SgmRqusmat->ptr[iq][jq]
          - LamLmat->ptr[iq][jq] - Grninvs->ptr[iq][jq];
        LamAmat->ptr[iq][jq] =-Hqus->ptr[iq][jq] + SgmLqusmat->ptr[iq][jq] - SgmAqusmat->ptr[iq][jq]
          + LamLmat->ptr[iq][jq] + Grninvs->ptr[Nqus+iq][Nqus+jq];
      }
      LamRmat->ptr[iq][iq] += Eng;
      LamAmat->ptr[iq][iq] += Eng;
    }
  }
  
  mat2d_del_dcmplx(Grninvs);
  mat2d_del_dcmplx(Hqus);
  
  return;
  
 error:
  abort();
}

static void dmft_pi(void)
{
  FILE  *fp = NULL;
  int   ispin, ig, sign;
  dreal dE, Eng, coef, BoseL, BoseG;
  
  dE = (Emax - Emin) / (Nmesh - 1);
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    for(ig = 0; ig < Nmesh; ++ig) {
      Eng  = Emin + dE*ig;
      sign = (Eng > 0.0) - (Eng < 0.0);
      if(sign != 0) {
        coef  = erf(fabs(Eng*Eng*Eng)/(OmegOptic*OmegOptic*OmegOptic)); 
        BoseL = coef * (1.0 / (exp(BetaOptic*fabs(Eng)) - 1.0) + 0.5 - 0.5*sign);
        BoseG = coef * (1.0 / (exp(BetaOptic*fabs(Eng)) - 1.0) + 0.5 + 0.5*sign);
      } else {
        coef  = 0.0;
        BoseL = coef;
        BoseG = coef;
      }
      PiL[ispin]->xval[ig] = Eng; PiG[ispin]->xval[ig] = Eng;
      PiR[ispin]->xval[ig] = Eng; PiA[ispin]->xval[ig] = Eng;
      PiL[ispin]->ptr[ig]->addr[0] =-I * RateOptic * BoseL;
      PiG[ispin]->ptr[ig]->addr[0] =-I * RateOptic * BoseG;
      PiR[ispin]->ptr[ig]->addr[0] =-0.5 * I * RateOptic * coef * sign;
      PiA[ispin]->ptr[ig]->addr[0] = 0.5 * I * RateOptic * coef * sign;
    }
  }
  
  if(Rank == Root) {
    fp = fopen("Pi.ascii", "w");
    for(ig = 0; ig < Nmesh; ++ig) {
      Eng = Emin + dE*ig;
      fprintf(fp, "%15.5e ", Eng);
      for(ispin = 0; ispin < Nspin; ++ispin)
        fprintf(fp, "%15.5e %15.5e %15.5e %15.5e ",
                cimag(PiR[ispin]->ptr[ig]->addr[0]),
                cimag(PiA[ispin]->ptr[ig]->addr[0]),
                cimag(PiL[ispin]->ptr[ig]->addr[0]),
                cimag(PiG[ispin]->ptr[ig]->addr[0]));
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

static void dmft_hybrid(void)
{
  int            ispin;
  matfunc_dcmplx *GrnRqusfunc = NULL, *GrnAqusfunc = NULL,
                 *GrnLqusfunc = NULL, *GrnGqusfunc = NULL,
                 *SgmRqusfunc = NULL, *SgmAqusfunc = NULL,
                 *SgmLqusfunc = NULL, *SgmGqusfunc = NULL,
                 *LamRfunc = NULL, *LamAfunc = NULL,
                 *LamLfunc = NULL, *LamGfunc = NULL,
                 *GrnMfunc = NULL;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnMfunc = GrnM[ispin];
    GrnRqusfunc = GrnRqus[ispin]; GrnAqusfunc = GrnAqus[ispin];
    GrnLqusfunc = GrnLqus[ispin]; GrnGqusfunc = GrnGqus[ispin];
    SgmRqusfunc = SgmRqus[ispin]; SgmAqusfunc = SgmAqus[ispin];
    SgmLqusfunc = SgmLqus[ispin]; SgmGqusfunc = SgmGqus[ispin];
    LamRfunc = LamR[ispin]; LamAfunc = LamA[ispin];
    LamLfunc = LamL[ispin]; LamGfunc = LamG[ispin];
    
    dmft_grnm(ispin, GrnMfunc);
    dmft_grnqus(GrnMfunc, GrnRqusfunc, GrnAqusfunc, GrnLqusfunc, GrnGqusfunc);
    dmft_lam(GrnRqusfunc, GrnAqusfunc, GrnLqusfunc, GrnGqusfunc,
             SgmRqusfunc, SgmAqusfunc, SgmLqusfunc, SgmGqusfunc,
             LamRfunc, LamAfunc, LamLfunc, LamGfunc);
  }
  
  if(isOptic) dmft_pi();
}

static void dmft_hybaug(void)
{
  int ispin;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    dmft_matfunc_aug(LamR[ispin], LamR_aug[ispin]);
    dmft_matfunc_aug(LamA[ispin], LamA_aug[ispin]);
    dmft_matfunc_aug(LamL[ispin], LamL_aug[ispin]);
    dmft_matfunc_aug(LamG[ispin], LamG_aug[ispin]);
    if(isOptic) {
      dmft_matfunc_aug(PiR[ispin], PiR_aug[ispin]);
      dmft_matfunc_aug(PiA[ispin], PiA_aug[ispin]);
      dmft_matfunc_aug(PiL[ispin], PiL_aug[ispin]);
      dmft_matfunc_aug(PiG[ispin], PiG_aug[ispin]);
    }
  }
}

static void dmft_sgm(const matfunc_dcmplx *LamRfunc, const matfunc_dcmplx *LamAfunc,
                     const matfunc_dcmplx *LamLfunc, const matfunc_dcmplx *LamGfunc,
                     const matfunc_dcmplx *GrnRqusfunc, const matfunc_dcmplx *GrnAqusfunc,
                     const matfunc_dcmplx *GrnLqusfunc, const matfunc_dcmplx *GrnGqusfunc,
                     matfunc_dcmplx *SgmRqusfunc, matfunc_dcmplx *SgmAqusfunc,
                     matfunc_dcmplx *SgmLqusfunc, matfunc_dcmplx *SgmGqusfunc)
{
  int          iq, jq, idim, jdim, ig;
  dreal        dE, Eng;
  mat2d_dcmplx *Hqus = NULL, *Grninvs = NULL;
  mat2d_dcmplx *GrnRqusmat = NULL, *GrnAqusmat = NULL,
               *GrnLqusmat = NULL, *GrnGqusmat = NULL,
               *SgmRqusmat = NULL, *SgmAqusmat = NULL,
               *SgmLqusmat = NULL, *SgmGqusmat = NULL,
               *LamRmat = NULL, *LamAmat = NULL,
               *LamLmat = NULL, *LamGmat = NULL;
  
  Hqus = mat2d_new_dcmplx();
  check_mem(Hqus, "Hqus");
  Hqus->alloc(Nqus, Nqus, Hqus);
  
  Grninvs = mat2d_new_dcmplx();
  check_mem(Grninvs, "Grninvs");
  Grninvs->alloc(NDIMSUB*Nqus, NDIMSUB*Nqus, Grninvs);
  
  for(iq = 0; iq < Nqus; ++iq) {
    idim = Indqus[iq];
    for(jq = 0; jq < Nqus; ++jq) {
      jdim = Indqus[jq];
      Hqus->ptr[iq][jq] = Hmat->ptr[idim][jdim];
    }
  }
  
  dE = (Emax - Emin) / (Nmesh - 1);
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    SgmRqusfunc->xval[ig] = Eng; SgmAqusfunc->xval[ig] = Eng;
    SgmLqusfunc->xval[ig] = Eng; SgmGqusfunc->xval[ig] = Eng;
    LamRmat = LamRfunc->ptr[ig]; LamAmat = LamAfunc->ptr[ig];
    LamLmat = LamLfunc->ptr[ig]; LamGmat = LamGfunc->ptr[ig];
    GrnRqusmat = GrnRqusfunc->ptr[ig]; GrnAqusmat = GrnAqusfunc->ptr[ig];
    GrnLqusmat = GrnLqusfunc->ptr[ig]; GrnGqusmat = GrnGqusfunc->ptr[ig];
    SgmRqusmat = SgmRqusfunc->ptr[ig]; SgmAqusmat = SgmAqusfunc->ptr[ig];
    SgmLqusmat = SgmLqusfunc->ptr[ig]; SgmGqusmat = SgmGqusfunc->ptr[ig];
    
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        Grninvs->ptr[iq][jq] = GrnLqusmat->ptr[iq][jq] + GrnRqusmat->ptr[iq][jq];
        Grninvs->ptr[Nqus+iq][jq] = GrnLqusmat->ptr[iq][jq];
        Grninvs->ptr[iq][Nqus+jq] = GrnGqusmat->ptr[iq][jq];
        Grninvs->ptr[Nqus+iq][Nqus+jq] = GrnLqusmat->ptr[iq][jq] - GrnAqusmat->ptr[iq][jq];
      }
    }
    
    lapack_zinvs(NDIMSUB*Nqus, Grninvs->addr);
    
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        SgmLqusmat->ptr[iq][jq] = Grninvs->ptr[Nqus+iq][jq] - LamLmat->ptr[iq][jq];
        SgmGqusmat->ptr[iq][jq] = Grninvs->ptr[iq][Nqus+jq] - LamGmat->ptr[iq][jq];
        SgmRqusmat->ptr[iq][jq] =-Hqus->ptr[iq][jq] - LamLmat->ptr[iq][jq] - LamRmat->ptr[iq][jq]
          - SgmLqusmat->ptr[iq][jq] - Grninvs->ptr[iq][jq];
        SgmAqusmat->ptr[iq][jq] =-Hqus->ptr[iq][jq] + LamLmat->ptr[iq][jq] - LamAmat->ptr[iq][jq]
          + SgmLqusmat->ptr[iq][jq] + Grninvs->ptr[Nqus+iq][Nqus+jq];
      }
      SgmRqusmat->ptr[iq][iq] += Eng;
      SgmAqusmat->ptr[iq][iq] += Eng;
    }
  }
  
  mat2d_del_dcmplx(Grninvs);
  mat2d_del_dcmplx(Hqus);
  
  return;
  
 error:
  abort();
}

static void dmft_selfen(void)
{
  int            ispin;
  matfunc_dcmplx *LamRfunc = NULL, *LamAfunc = NULL,
                 *LamLfunc = NULL, *LamGfunc = NULL,
                 *GrnRqusfunc = NULL, *GrnAqusfunc = NULL,
                 *GrnLqusfunc = NULL, *GrnGqusfunc = NULL,
                 *SgmRqusfunc = NULL, *SgmAqusfunc = NULL,
                 *SgmLqusfunc = NULL, *SgmGqusfunc = NULL;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    LamRfunc = LamR[ispin]; LamAfunc = LamA[ispin];
    LamLfunc = LamL[ispin]; LamGfunc = LamG[ispin];
    GrnRqusfunc = GrnRqus[ispin]; GrnAqusfunc = GrnAqus[ispin];
    GrnLqusfunc = GrnLqus[ispin]; GrnGqusfunc = GrnGqus[ispin];
    SgmRqusfunc = SgmRqus[ispin]; SgmAqusfunc = SgmAqus[ispin];
    SgmLqusfunc = SgmLqus[ispin]; SgmGqusfunc = SgmGqus[ispin];
    
    dmft_sgm(LamRfunc, LamAfunc, LamLfunc, LamGfunc,
             GrnRqusfunc, GrnAqusfunc, GrnLqusfunc, GrnGqusfunc,
             SgmRqusfunc, SgmAqusfunc, SgmLqusfunc, SgmGqusfunc);
  }
}

void dmft_sumdosqus(matfunc_dcmplx * const *GrnRqusfunc, dreal **SumDosqus)
{
  int   ispin, ig, iq;
  dreal *xval = NULL, *dosqus = NULL;
  
  dosqus = (dreal *) calloc(Nmesh, sizeof(dreal));
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    xval = GrnRqusfunc[ispin]->xval;
    for(iq = 0; iq < Nqus; ++iq) {
      for(ig = 0; ig < Nmesh; ++ig) {
        dosqus[ig] =-cimag(GrnRqusfunc[ispin]->ptr[ig]->ptr[iq][iq]) / PI;
      }
      SumDosqus[ispin][iq] = numeric_quadrat_dreal(Nmesh, xval, dosqus);
    }
  }
  
  freeup(dosqus);
}

void dmft_sumpopqus(matfunc_dcmplx * const *GrnLqusfunc, dreal **SumPopqus)
{
  int   ispin, ig, iq;
  dreal *xval = NULL, *popqus = NULL;
  
  popqus = (dreal *) calloc(Nmesh, sizeof(dreal));
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    xval = GrnLqusfunc[ispin]->xval;
    for(iq = 0; iq < Nqus; ++iq) {
      for(ig = 0; ig < Nmesh; ++ig) {
        popqus[ig] = cimag(GrnLqusfunc[ispin]->ptr[ig]->ptr[iq][iq]) / (2.0*PI);
      }
      SumPopqus[ispin][iq] = numeric_quadrat_dreal(Nmesh, xval, popqus);
    }
  }
  
  freeup(popqus);
}

static int dmft_cnvg(void)
{
  int   ispin, iq, cnvg;
  dreal popdiff, maxpopdiff;
  
  cnvg = 0;
  
  for(ispin = 0, maxpopdiff = 0.0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      popdiff = fabs(SumPopqus_new[ispin][iq] - SumPopqus_old[ispin][iq]);
      if(popdiff > maxpopdiff) maxpopdiff = popdiff;
    }
  }
  
  if(Rank == Root) {
    printf("%s\n", dmftlabel);
    printf("%s maxpopdiff = %15.5e\n", dmftlabel, maxpopdiff);
    printf("%s\n", dmftlabel);
  }
  
  if(maxpopdiff < TOLDPOP) cnvg = 1;
  
  return cnvg;
}

static void dmft_dosqus_out(void)
{
  FILE  *fp = NULL;
  int   ispin, iq, ngrd, ig;
  dreal **SumDosqus = NULL;
  char  *filename = NULL;
  
  SumDosqus = (dreal **) calloc(Nspin, sizeof(dreal *));
  for(ispin = 0; ispin < Nspin; ++ispin)
    SumDosqus[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
  
  dmft_sumdosqus(GrnRqus, SumDosqus);
  
  if(Rank == Root) {
    printf("%s\n", dmftlabel);
    printf("%s SumDosqus =\n", dmftlabel);
    for(ispin = 0; ispin < Nspin; ++ispin)
      for(iq = 0; iq < Nqus; ++iq)
        printf("%s %5d %15.5e\n", dmftlabel, ispin, SumDosqus[ispin][iq]);
  }
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      ngrd = GrnRqus[ispin]->ngrd;
      if(Dosqus[ispin][iq]) freeup(Dosqus[ispin][iq]);
      Dosqus[ispin][iq] = (dreal *) calloc(ngrd, sizeof(dreal));
      check_mem(Dosqus[ispin][iq], "Dosqus[ispin][iq]");
      for(ig = 0; ig < ngrd; ++ig) {
        Dosqus[ispin][iq][ig] =-cimag(GrnRqus[ispin]->ptr[ig]->ptr[iq][iq]) / PI;
      }
      // Output
      filename = (char *) calloc(SHRT_MAX, sizeof(char));
      sprintf(filename, "%s_%1d_%04d.ascii", "Dosqus", ispin, iq);
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", GrnRqus[ispin]->xval[ig], Dosqus[ispin][iq][ig]);
        }
        fclose(fp);
      }
      freeup(filename);
    }  
  }
  
  for(ispin = Nspin - 1; ispin > -1; --ispin)
    freeup(SumDosqus[ispin]);
  freeup(SumDosqus);
  
  return;
  
 error:
  abort();
}

static void dmft_popqus_out(void)
{
  FILE  *fp = NULL;
  int   ispin, iq, ngrd, ig;
  dreal **SumPopqus = NULL;
  char  *filename = NULL;
  
  SumPopqus = (dreal **) calloc(Nspin, sizeof(dreal *));
  for(ispin = 0; ispin < Nspin; ++ispin)
    SumPopqus[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
  
  dmft_sumpopqus(GrnLqus, SumPopqus);
  
  if(Rank == Root) {
    printf("%s\n", dmftlabel);
    printf("%s SumPopqus =\n", dmftlabel);
    for(ispin = 0; ispin < Nspin; ++ispin)
      for(iq = 0; iq < Nqus; ++iq)
        printf("%s %5d %15.5e\n", dmftlabel, ispin, SumPopqus[ispin][iq]);
  }
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      ngrd = GrnLqus[ispin]->ngrd;
      if(Popqus[ispin][iq]) freeup(Popqus[ispin][iq]);
      Popqus[ispin][iq] = (dreal *) calloc(ngrd, sizeof(dreal));
      check_mem(Popqus[ispin][iq], "Popqus[ispin][iq]");
      for(ig = 0; ig < ngrd; ++ig) {
        Popqus[ispin][iq][ig] = cimag(GrnLqus[ispin]->ptr[ig]->ptr[iq][iq]) / (2.0*PI);
      }
      // Output
      filename = (char *) calloc(SHRT_MAX, sizeof(char));
      sprintf(filename, "%s_%1d_%04d.ascii", "Popqus", ispin, iq);
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", GrnLqus[ispin]->xval[ig], Popqus[ispin][iq][ig]);
        }
        fclose(fp);
      }
      freeup(filename);
    }  
  }
  
  for(ispin = Nspin - 1; ispin > -1; --ispin)
    freeup(SumPopqus[ispin]);
  freeup(SumPopqus);
  
  return;
  
 error:
  abort();
}

static void dmft_itot(void)
{
  int          ienv, ispin, idim, idimsub, ig;
  dreal        dE, Eng;
  dcmplx       grnL, grnG, senL, senG;
  dreal        *xval = NULL;
  dcmplx       *cfun = NULL;
  mat2d_dcmplx *GrnMmat = NULL;
  
  xval = (dreal *) calloc(Nmesh, sizeof(dreal));
  cfun = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  
  dE = (Emax - Emin) / (Nmesh - 1);
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin + dE * ig;
    xval[ig] = Eng;
  }
  
  for(ienv = 0, Itot[ienv] = 0.0; ienv < Nenv; ++ienv) {
    for(ispin = 0; ispin < Nspin; ++ispin) {
      for(idim = 0; idim < Ndim; ++idim) {
        idimsub = NDIMSUB * idim;
        for(ig = 0; ig < Nmesh; ++ig) {
          Eng = xval[ig];
          GrnMmat= GrnM[ispin]->ptr[ig];
          grnL = GrnMmat->ptr[idimsub+1][idimsub];
          grnG = GrnMmat->ptr[idimsub][idimsub+1];
          senL = I * Gamma[ienv][idim] * stats_ferm(Beta, Muenv[ienv], Eng);
          senG =-I * Gamma[ienv][idim] * stats_invsferm(Beta, Muenv[ienv], Eng);
          cfun[ig] = grnG * senL - grnL * senG;
        }
        Itot[ienv] += numeric_quadrat_dcmplx(Nmesh, xval, cfun) / (2.0*PI);
      }
    }
  }
  
  if(Rank == Root) {
    printf("%s\n", dmftlabel);
    printf("%s ienv, Itot =\n", dmftlabel);
    for(ienv = 0; ienv < Nenv; ++ienv) {
      printf("%s %5d %15.5e %15.5e\n", dmftlabel, ienv,
             creal(Itot[ienv]), cimag(Itot[ienv]));
    }
  }
  
  freeup(cfun); freeup(xval);
}

static void dmft_mixing_grnqus(void)
{
  int            ispin, ig, iq, jq;
  mat2d_dcmplx   *GrnRqusmat = NULL, *GrnAqusmat = NULL, *GrnLqusmat = NULL, *GrnGqusmat = NULL;
  mat2d_dcmplx   *GrnRqusmatPT = NULL, *GrnAqusmatPT = NULL, *GrnLqusmatPT = NULL, *GrnGqusmatPT = NULL;
  matfunc_dcmplx *GrnRqusfunc = NULL, *GrnAqusfunc = NULL, *GrnLqusfunc = NULL, *GrnGqusfunc = NULL;
  matfunc_dcmplx *GrnRqusfuncPT = NULL, *GrnAqusfuncPT = NULL, *GrnLqusfuncPT = NULL, *GrnGqusfuncPT = NULL;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnRqusfunc = GrnRqus[ispin]; GrnAqusfunc = GrnAqus[ispin];
    GrnLqusfunc = GrnLqus[ispin]; GrnGqusfunc = GrnGqus[ispin];
    
    GrnRqusfuncPT = GrnRqusPT[ispin]; GrnAqusfuncPT = GrnAqusPT[ispin];
    GrnLqusfuncPT = GrnLqusPT[ispin]; GrnGqusfuncPT = GrnGqusPT[ispin];
    
    for(ig = 0; ig < Nmesh; ++ig) {
      GrnRqusmat = GrnRqusfunc->ptr[ig]; GrnAqusmat = GrnAqusfunc->ptr[ig];
      GrnLqusmat = GrnLqusfunc->ptr[ig]; GrnGqusmat = GrnGqusfunc->ptr[ig];
      
      GrnRqusmatPT = GrnRqusfuncPT->ptr[ig]; GrnAqusmatPT = GrnAqusfuncPT->ptr[ig];
      GrnLqusmatPT = GrnLqusfuncPT->ptr[ig]; GrnGqusmatPT = GrnGqusfuncPT->ptr[ig];
      for(iq = 0; iq < Nqus; ++iq) {
        for(jq = 0; jq < Nqus; ++jq) {
          GrnRqusmat->ptr[iq][jq] = Mixing * GrnRqusmat->ptr[iq][jq] + (1.0-Mixing) * GrnRqusmatPT->ptr[iq][jq];
          GrnAqusmat->ptr[iq][jq] = Mixing * GrnAqusmat->ptr[iq][jq] + (1.0-Mixing) * GrnAqusmatPT->ptr[iq][jq];
          GrnLqusmat->ptr[iq][jq] = Mixing * GrnLqusmat->ptr[iq][jq] + (1.0-Mixing) * GrnLqusmatPT->ptr[iq][jq];
          GrnGqusmat->ptr[iq][jq] = Mixing * GrnGqusmat->ptr[iq][jq] + (1.0-Mixing) * GrnGqusmatPT->ptr[iq][jq];
        }
      }
    }
  }
}

void dmft_main(void)
{
  int cnvg, iter, ispin, iq;
  
  if(Rank == Root) {
    printf("%s\n", dmftlabel);
    printf("%s *********** DMFT Loop **********\n", dmftlabel);
    fflush(stdout);
  }
  
  for(cnvg = 0, iter = 0; iter < MAXITER && !cnvg; ++iter) {
    if(Rank == Root) {
      printf("%s ITER %5d\n", dmftlabel, iter); 
      fflush(stdout);
    }
    // Hybridization for impurity
    dmft_hybrid(); dmft_hybaug();
    dmft_sumpopqus(GrnLqus, SumPopqus_old);
    if(Rank == Root) {
      printf("%s SumPopqus_old =\n", dmftlabel);
      for(ispin = 0; ispin < Nspin; ++ispin)
        for(iq = 0; iq < Nqus; ++iq)
          printf("%s %5d %5d %15.5e\n", dmftlabel, ispin, iq, SumPopqus_old[ispin][iq]);
    }
    // NCA impurity solver
    if(Mixing > THRSHLD) {
      nca_psp_rscf();
      nca_psp_gprj();
      nca_psp_aprj();
      nca_psp_lscf();
      nca_psp_out();
      nca_grnqus();
    }
    // End of NCA
    // Perturbation impurity solver
    if(1.0-Mixing > THRSHLD) perturb_scf();
    // Mixing of NCA and perturbation
    dmft_mixing_grnqus();
    //
    dmft_sumpopqus(GrnLqus, SumPopqus_new);
    if(Rank == Root) {
      printf("%s SumPopqus_new =\n", dmftlabel);
      for(ispin = 0; ispin < Nspin; ++ispin)
        for(iq = 0; iq < Nqus; ++iq)
          printf("%s %5d %5d %15.5e\n", dmftlabel, ispin, iq, SumPopqus_new[ispin][iq]);
    }
    // Convergence check
    cnvg = dmft_cnvg();
    // Update self energies
    dmft_selfen();
    
    if(Rank == Root) {
      printf("%s ********************************\n", dmftlabel);
      fflush(stdout);
    }
  }
  
  dmft_itot();
  dmft_dosqus_out();
  dmft_popqus_out();
}
