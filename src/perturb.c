#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "dbg.h"
#include "sys.h"
#include "lapack.h"
#include "numeric.h"
#include "constants.h"
#include "variables.h"
#include "dmft.h"

static const int      MAXNSCF = 500;
static const dreal    MIX     = 9e-1;
static const dreal    TOLDPOP = 1e-9;
static dreal          *Vhar = NULL, *UU = NULL;
static dreal          **Occqus_old = NULL, **Occqus_new = NULL;
static mat2d_dcmplx   *Hpert = NULL;
static matfunc_dcmplx **GrnMqus = NULL;

static const char *perturbin = "perturb.in";
static const char *ptlabel   = "perturb:        ";

void perturb_alloc(void)
{
  FILE  *fp = NULL;
  int   iq, jq, ispin;
  dreal tmpre, tmpim;
  char  line[SHRT_MAX];
  // Hpert
  Hpert = mat2d_new_dcmplx();
  check_mem(Hpert, "Hpert");
  Hpert->alloc(Nqus, Nqus, Hpert);
  // Hartree term
  Vhar = (dreal *) calloc(Nqus, sizeof(dreal));
  // Local energy and Hubbard parameter
  UU = (dreal *) calloc(Nqus, sizeof(dreal));
  
  if(Rank == Root) {
    fp = fopen(perturbin, "r");
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf %lf", &tmpre, &tmpim);
        Hpert->ptr[iq][jq] = tmpre + I * tmpim;
      }
    }
    for(iq = 0; iq < Nqus; ++iq)
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", UU+iq);
    fclose(fp);
  }
  
  sys_bcast_dreal(UU, Nqus, Root);
  sys_bcast_dcmplx(Hpert->addr, Nqus*Nqus, Root);
  
  for(iq = 0; iq < Nqus; ++iq) {
    for(jq = 0; jq < Nqus; ++jq)
      Hpert->ptr[iq][jq] /= EV_HA;
    UU[iq] /= EV_HA;
  }
  
  // Population
  Occqus_old = (dreal **) calloc(Nspin, sizeof(dreal *));
  Occqus_new = (dreal **) calloc(Nspin, sizeof(dreal *));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    Occqus_old[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
    Occqus_new[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
  }
  // Perturbation Green functions
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnRqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(GrnRqusPT[ispin], "GrnRqusPT[ispin]");
    GrnRqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, GrnRqusPT[ispin]);
    
    GrnAqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(GrnAqusPT[ispin], "GrnAqusPT[ispin]");
    GrnAqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, GrnAqusPT[ispin]);
    
    GrnLqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(GrnLqusPT[ispin], "GrnLqusPT[ispin]");
    GrnLqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, GrnLqusPT[ispin]);
    
    GrnGqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(GrnGqusPT[ispin], "GrnGqusPT[ispin]");
    GrnGqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, GrnGqusPT[ispin]);
  }
  // Perturbation self energies
  for(ispin = 0; ispin < Nspin; ++ispin) {
    SgmRqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(SgmRqusPT[ispin], "SgmRqusPT[ispin]");
    SgmRqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, SgmRqusPT[ispin]);
    
    SgmAqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(SgmAqusPT[ispin], "SgmAqusPT[ispin]");
    SgmAqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, SgmAqusPT[ispin]);
    
    SgmLqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(SgmLqusPT[ispin], "SgmLqusPT[ispin]");
    SgmLqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, SgmLqusPT[ispin]);
    
    SgmGqusPT[ispin] = matfunc_new_dcmplx();
    check_mem(SgmGqusPT[ispin], "SgmGqusPT[ispin]");
    SgmGqusPT[ispin]->alloc(Nmesh, Nqus, Nqus, SgmGqusPT[ispin]);
  }
  
  GrnMqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnMqus[ispin] = matfunc_new_dcmplx();
    check_mem(GrnMqus[ispin], "GrnMqus[ispin]");
    GrnMqus[ispin]->alloc(Nmesh, NDIMSUB*Nqus, NDIMSUB*Nqus, GrnMqus[ispin]);
  }
  
  return;
  
 error:
  abort();
}

void perturb_deall(void)
{
  int ispin;
  
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    matfunc_del_dcmplx(GrnMqus[ispin]);
  }
  freeup(GrnMqus);
  
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    matfunc_del_dcmplx(SgmGqusPT[ispin]);
    matfunc_del_dcmplx(SgmLqusPT[ispin]);
    matfunc_del_dcmplx(SgmAqusPT[ispin]);
    matfunc_del_dcmplx(SgmRqusPT[ispin]);
  }
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    matfunc_del_dcmplx(GrnGqusPT[ispin]);
    matfunc_del_dcmplx(GrnLqusPT[ispin]);
    matfunc_del_dcmplx(GrnAqusPT[ispin]);
    matfunc_del_dcmplx(GrnRqusPT[ispin]);
  }
  
  for(ispin = Nspin - 1; ispin > -1; --ispin) {
    freeup(Occqus_new[ispin]);
    freeup(Occqus_old[ispin]);
  }
  freeup(Occqus_new);
  freeup(Occqus_old);
  
  mat2d_del_dcmplx(Hpert);
  
  freeup(UU);
  freeup(Vhar);
}

void perturb_print(void)
{
  int iq;
  
  if(Rank == Root) {
    printf("%s\n", ptlabel);
    printf("%s UU (Ha) =\n", ptlabel);
    for(iq = 0; iq < Nqus; ++iq)
      printf("%s %5d %15.5e\n", ptlabel, iq, UU[iq]);
    printf("%s\n", ptlabel);
    fflush(stdout);
  }
}

static void perturb_copy_grnqus(void)
{
  int ispin;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    matfunc_copy_dcmplx(GrnRqus[ispin], GrnRqusPT[ispin]);
    matfunc_copy_dcmplx(GrnAqus[ispin], GrnAqusPT[ispin]);
    matfunc_copy_dcmplx(GrnLqus[ispin], GrnLqusPT[ispin]);
    matfunc_copy_dcmplx(GrnGqus[ispin], GrnGqusPT[ispin]);
  }
}

static void perturb_grnmqus(int ispin, matfunc_dcmplx *GrnMfunc)
{
  int            ig, iq, jq, iqsub, jqsub;
  dreal          dE, Eng;
  mat2d_dcmplx   *GrnMmat = NULL;
  mat2d_dcmplx   *LamRmat = NULL, *LamAmat = NULL, *LamLmat = NULL, *LamGmat = NULL,
                 *SgmRmat = NULL, *SgmAmat = NULL, *SgmLmat = NULL, *SgmGmat = NULL;
  matfunc_dcmplx *LamRfunc = NULL, *LamAfunc = NULL, *LamLfunc = NULL, *LamGfunc = NULL,
                 *SgmRfunc = NULL, *SgmAfunc = NULL, *SgmLfunc = NULL, *SgmGfunc = NULL;
  
  dE = (Emax - Emin) / (Nmesh - 1);
  
  LamRfunc = LamR[ispin]; LamAfunc = LamA[ispin];
  LamLfunc = LamL[ispin]; LamGfunc = LamG[ispin];
  
  SgmRfunc = SgmRqusPT[ispin]; SgmAfunc = SgmAqusPT[ispin];
  SgmLfunc = SgmLqusPT[ispin]; SgmGfunc = SgmGqusPT[ispin];
  
  for(ig = 0; ig < Nmesh; ++ig) {
    Eng = Emin+ dE * ig;
    GrnMfunc->xval[ig] = Eng;
    GrnMmat = GrnMfunc->ptr[ig];
    
    LamRmat = LamRfunc->ptr[ig]; LamAmat = LamAfunc->ptr[ig];
    LamLmat = LamLfunc->ptr[ig]; LamGmat = LamGfunc->ptr[ig];
    
    SgmRmat = SgmRfunc->ptr[ig]; SgmAmat = SgmAfunc->ptr[ig];
    SgmLmat = SgmLfunc->ptr[ig]; SgmGmat = SgmGfunc->ptr[ig];
    // Loop over subblock
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq) {
        // Part 1
        jqsub = NDIMSUB * jq;
        iqsub = NDIMSUB * iq;
        GrnMmat->ptr[iqsub][jqsub] =-Hpert->ptr[iq][jq] - LamLmat->ptr[iq][jq] - LamRmat->ptr[iq][jq]
          - SgmLmat->ptr[iq][jq] - SgmRmat->ptr[iq][jq];
        // Part 2
        jqsub = NDIMSUB * jq;
        iqsub = NDIMSUB * iq + 1;
        GrnMmat->ptr[iqsub][jqsub] = LamLmat->ptr[iq][jq] + SgmLmat->ptr[iq][jq];
        // Part 3
        jqsub = NDIMSUB * jq + 1;
        iqsub = NDIMSUB * iq;
        GrnMmat->ptr[iqsub][jqsub] = LamGmat->ptr[iq][jq] + SgmGmat->ptr[iq][jq];
        // Part 4
        jqsub = NDIMSUB * jq + 1;
        iqsub = NDIMSUB * iq + 1;
        GrnMmat->ptr[iqsub][jqsub] = Hpert->ptr[iq][jq] - LamLmat->ptr[iq][jq] + LamAmat->ptr[iq][jq]
          - SgmLmat->ptr[iq][jq] + SgmAmat->ptr[iq][jq];
      }
      // Diagonal part
      iqsub = NDIMSUB * iq;
      GrnMmat->ptr[iqsub][iqsub] += Eng - Vhar[iq];
      iqsub = NDIMSUB * iq + 1;
      GrnMmat->ptr[iqsub][iqsub] +=-Eng + Vhar[iq];
    }
    lapack_zinvs(NDIMSUB*Nqus, GrnMmat->addr);
  }
}

static void perturb_vhar(void)
{
  int   ispin, iq;
  dreal **Occ = NULL;
  
  Occ = (dreal **) calloc(Nspin, sizeof(dreal *));
  for(ispin = 0; ispin < Nspin; ++ispin)
    Occ[ispin] = (dreal *) calloc(Nqus, sizeof(dreal));
  
  dmft_sumpopqus(GrnLqusPT, Occ);
  
  for(iq = 0; iq < Nqus; ++iq)
    Vhar[iq] = UU[iq] * Occ[0][iq];
  
  for(ispin = Nspin - 1; ispin > -1; --ispin)
    freeup(Occ[ispin]);
  freeup(Occ);
}

static void perturb_sgmqus(matfunc_dcmplx *GrnRqusfunc, matfunc_dcmplx *GrnAqusfunc,
                           matfunc_dcmplx *GrnLqusfunc, matfunc_dcmplx *GrnGqusfunc,
                           matfunc_dcmplx *SgmRqusfunc, matfunc_dcmplx *SgmAqusfunc,
                           matfunc_dcmplx *SgmLqusfunc, matfunc_dcmplx *SgmGqusfunc)
{
  int    iq, ig;
  dreal  dE, dT;
  dcmplx *grnLf = NULL, *grnGf = NULL, *grnLt = NULL, *grnGt = NULL,
         *sgmLf = NULL, *sgmGf = NULL, *sgmLt = NULL, *sgmGt = NULL,
         *fa = NULL, *fb = NULL;
  
  dE = (Emax - Emin) / (Nmesh - 1);
  dT = 2.0*PI / (Emax - Emin);
  /*
  for(ig = 0; ig < Nmesh; ++ig) {
    SgmRqusfunc->ptr[ig]->reset(SgmRqusfunc->ptr[ig]);
    SgmAqusfunc->ptr[ig]->reset(SgmAqusfunc->ptr[ig]);
    SgmLqusfunc->ptr[ig]->reset(SgmLqusfunc->ptr[ig]);
    SgmGqusfunc->ptr[ig]->reset(SgmGqusfunc->ptr[ig]);
  }
  */
  grnLf = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  grnGf = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  grnLt = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  grnGt = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  
  sgmLf = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  sgmGf = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  sgmLt = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  sgmGt = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  
  fa = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  fb = (dcmplx *) calloc(Nmesh, sizeof(dcmplx));
  
  for(iq = 0; iq < Nqus; ++iq) {
    for(ig = 0; ig < Nmesh; ++ig) {
      grnLf[ig] = GrnLqusfunc->ptr[ig]->ptr[iq][iq];
      grnGf[ig] = GrnGqusfunc->ptr[ig]->ptr[iq][iq];
    }
    numeric_freq2time(Nmesh-1, dE, grnLf, grnLt); grnLt[Nmesh-1] = grnLt[0];
    numeric_freq2time(Nmesh-1, dE, grnGf, grnGt); grnGt[Nmesh-1] = grnGt[0];
    // Lesser
    numeric_mirror(Nmesh-1, grnGt, fa); fa[Nmesh-1] = fa[0];
    for(ig = 0; ig < Nmesh; ++ig)
      sgmLt[ig] = UU[iq]*UU[iq]*grnLt[ig]*grnLt[ig]*fa[ig];
    numeric_time2freq(Nmesh-1, dT, sgmLt, sgmLf); sgmLf[Nmesh-1] = sgmLf[0];
    for(ig = 0; ig < Nmesh; ++ig)
      SgmLqusfunc->ptr[ig]->ptr[iq][iq] = MIX*SgmLqusfunc->ptr[ig]->ptr[iq][iq] + (1.0-MIX) * sgmLf[ig];
    // Greater
    numeric_mirror(Nmesh-1, grnLt, fa); fa[Nmesh-1] = fa[0];
    for(ig = 0; ig < Nmesh; ++ig)
      sgmGt[ig] = UU[iq]*UU[iq]*grnGt[ig]*grnGt[ig]*fa[ig];
    numeric_time2freq(Nmesh-1, dT, sgmGt, sgmGf); sgmGf[Nmesh-1] = sgmGf[0];
    for(ig = 0; ig < Nmesh; ++ig)
      SgmGqusfunc->ptr[ig]->ptr[iq][iq] = MIX*SgmGqusfunc->ptr[ig]->ptr[iq][iq] + (1.0-MIX) * sgmGf[ig];
    // Retarded
    for(ig = 0; ig < Nmesh; ++ig)
      fa[ig] = sgmGf[ig] - sgmLf[ig];
    numeric_kronig(Nmesh-1, fa, fb);
    for(ig = 0; ig < Nmesh; ++ig) {
      SgmRqusfunc->ptr[ig]->ptr[iq][iq] = MIX*SgmRqusfunc->ptr[ig]->ptr[iq][iq] + (1.0-MIX) * fb[ig];
      SgmAqusfunc->ptr[ig]->ptr[iq][iq] = MIX*SgmAqusfunc->ptr[ig]->ptr[iq][iq] + (1.0-MIX) * conj(fb[ig]);
    }
  }
  
  freeup(fb); freeup(fa);
  
  freeup(sgmGt); freeup(sgmLt);
  freeup(sgmGf); freeup(sgmLf);
  freeup(grnGt); freeup(grnLt);
  freeup(grnGf); freeup(grnLf);
}

static void perturb_grnqus(matfunc_dcmplx *GrnMfunc, 
                           matfunc_dcmplx *GrnRqusfunc, matfunc_dcmplx *GrnAqusfunc,
                           matfunc_dcmplx *GrnLqusfunc, matfunc_dcmplx *GrnGqusfunc)
{
  int          ig, iq, jq, iqsub, jqsub;
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
      iqsub = NDIMSUB * iq;
      for(jq = 0; jq < Nqus; ++jq) {
        jqsub = NDIMSUB * jq;
        GrnLqusmat->ptr[iq][jq] = GrnMmat->ptr[iqsub+1][jqsub];
        GrnGqusmat->ptr[iq][jq] = GrnMmat->ptr[iqsub][jqsub+1];
        GrnRqusmat->ptr[iq][jq] = GrnMmat->ptr[iqsub][jqsub] - GrnMmat->ptr[iqsub+1][jqsub];
        GrnAqusmat->ptr[iq][jq] = GrnMmat->ptr[iqsub+1][jqsub] - GrnMmat->ptr[iqsub+1][jqsub+1];
      }
    }
  }
}

static void perturb_main(void)
{
  int            ispin;
  matfunc_dcmplx *GrnRfunc = NULL, *GrnAfunc = NULL, *GrnLfunc = NULL, *GrnGfunc = NULL,
                 *SgmRfunc = NULL, *SgmAfunc = NULL, *SgmLfunc = NULL, *SgmGfunc = NULL,
                 *GrnMfunc = NULL;
  
  perturb_vhar();
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    GrnMfunc = GrnMqus[ispin];
    
    GrnRfunc = GrnRqusPT[ispin]; GrnAfunc = GrnAqusPT[ispin];
    GrnLfunc = GrnLqusPT[ispin]; GrnGfunc = GrnGqusPT[ispin];
    
    SgmRfunc = SgmRqusPT[ispin]; SgmAfunc = SgmAqusPT[ispin];
    SgmLfunc = SgmLqusPT[ispin]; SgmGfunc = SgmGqusPT[ispin];
    
    perturb_sgmqus(GrnRfunc, GrnAfunc, GrnLfunc, GrnGfunc, SgmRfunc, SgmAfunc, SgmLfunc, SgmGfunc);
    perturb_grnmqus(ispin, GrnMfunc);
    perturb_grnqus(GrnMfunc, GrnRfunc, GrnAfunc, GrnLfunc, GrnGfunc);
  }
}

static void perturb_dosqus_out(void)
{
  FILE  *fp = NULL;
  int   ispin, iq, ngrd, ig;
  char  *filename = NULL;
  dreal *dos = NULL;
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      ngrd = GrnRqusPT[ispin]->ngrd;
      dos = (dreal *) calloc(ngrd, sizeof(dreal));
      for(ig = 0; ig < ngrd; ++ig) {
        dos[ig] =-cimag(GrnRqusPT[ispin]->ptr[ig]->ptr[iq][iq]) / PI;
      }
      // Output
      filename = (char *) calloc(SHRT_MAX, sizeof(char));
      sprintf(filename, "%s_%1d_%04d.ascii", "DosqusPT", ispin, iq);
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", GrnRqusPT[ispin]->xval[ig], dos[ig]);
        }
        fclose(fp);
      }
      freeup(filename);
      freeup(dos);
    }
  }
}

static void perturb_popqus_out(void)
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
        pop[ig] = cimag(GrnLqusPT[ispin]->ptr[ig]->ptr[iq][iq]) / (2.0*PI);
      }
      // Output
      filename = (char *) calloc(SHRT_MAX, sizeof(char));
      sprintf(filename, "%s_%1d_%04d.ascii", "PopqusPT", ispin, iq);
      if(Rank == Root) {
        fp = fopen(filename, "w");
        for(ig = 0; ig < ngrd; ++ig) {
          fprintf(fp, "%15.5e %15.5e\n", GrnLqusPT[ispin]->xval[ig], pop[ig]);
        }
        fclose(fp);
      }
      freeup(filename);
      freeup(pop);
    }  
  }
} 

static int perturb_cnvg(void)
{
  int   ispin, iq, cnvg;
  dreal popdiff, maxpopdiff;
  
  cnvg = 0;
  
  for(ispin = 0, maxpopdiff = 0.0; ispin < Nspin; ++ispin) {
    for(iq = 0; iq < Nqus; ++iq) {
      popdiff = fabs(Occqus_new[ispin][iq] - Occqus_old[ispin][iq]);
      if(popdiff > maxpopdiff) maxpopdiff = popdiff;
    }
  }
  
  if(Rank == Root) {
    printf("%s\n", ptlabel);
    printf("%s maxpopdiff = %15.5e\n", ptlabel, maxpopdiff);
    printf("%s\n", ptlabel);
  }
  
  if(maxpopdiff < TOLDPOP) cnvg = 1;
  
  return cnvg;
}

void perturb_scf(void)
{
  FILE *fp = NULL;
  int iscf, cnvg, ispin, iq;
  int ig, jq;
  
  perturb_copy_grnqus();
  
  fp = fopen("Lam.ascii", "w");
  for(ig = 0; ig < Nmesh; ++ig) {
    fprintf(fp, "%15.5e ", LamR[0]->xval[ig]);
    
    for(iq = 0; iq < Nqus; ++iq) {
      for(jq = 0; jq < Nqus; ++jq)
        fprintf(fp, "%15.5e %15.5e ", creal(LamR[0]->ptr[ig]->ptr[iq][jq]), cimag(LamR[0]->ptr[ig]->ptr[iq][jq]));
    }
    
    fprintf(fp, "\n");
  }
  fclose(fp);
  
  if(Rank == Root) printf("%s ========= Perturbation =========\n", ptlabel);
  for(iscf = 0, cnvg = 0; iscf < MAXNSCF && !cnvg; ++iscf) {
    if(Rank == Root) printf("%s ISCF %5d\n", ptlabel, iscf);
    dmft_sumpopqus(GrnLqusPT, Occqus_old);
    
    if(Rank == Root) {
      printf("%s Occqus_old =\n", ptlabel);
      for(ispin = 0; ispin < Nspin; ++ispin)
        for(iq = 0; iq < Nqus; ++iq)
          printf("%s %5d %5d %15.5e\n", ptlabel, ispin, iq, Occqus_old[ispin][iq]);
    }
    
    perturb_main();
    dmft_sumpopqus(GrnLqusPT, Occqus_new);
    if(Rank == Root) {
      printf("%s Occqus_new =\n", ptlabel);
      for(ispin = 0; ispin < Nspin; ++ispin)
        for(iq = 0; iq < Nqus; ++iq)
          printf("%s %5d %5d %15.5e\n", ptlabel, ispin, iq, Occqus_new[ispin][iq]);
    }
    cnvg = perturb_cnvg();
    if(Rank == Root) printf("%s --------------------------------\n", ptlabel);
  }
  if(Rank == Root) printf("%s ================================\n", ptlabel);
  
  perturb_dosqus_out();
  perturb_popqus_out();
}
