#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "dbg.h"
#include "sys.h"
#include "constants.h"
#include "numeric.h"
#include "lapack.h"
#include "variables.h"

int   Nmesh, Nenv, Nspin, Ndim, Nqus, Npsp, Nblk, Nmesh_aug;
int   isOptic, Gstat, Xstat;
dreal Efrm, Temp, Beta, Emin, Emax, Emin_aug, Emax_aug, Mixing;
dreal TempOptic, BetaOptic, OmegOptic, RateOptic;

int            *Npspblk = NULL, *Zeta = NULL, *Blkpsp = NULL, *Indqus = NULL;
int            **Pspblk = NULL;
dreal          *Muenv = NULL, *Epsp = NULL,
               *Sumdospsp = NULL, *Sumpoppsp = NULL;
dreal          **Gamma = NULL, **Dospsp = NULL, **Poppsp = NULL;
dreal          ***Dosqus = NULL, ***Popqus = NULL;
dcmplx         *Itot = NULL;
mat2d_dcmplx   *Hmat = NULL;
mat2d_dcmplx   **Hpsp = NULL;
mat2d_dcmplx   ***Xipsp = NULL;
mat2d_dcmplx   **Chipsp = NULL;
matfunc_dcmplx **PiR = NULL, **PiA = NULL, **PiL = NULL, **PiG = NULL;
matfunc_dcmplx **SgmRpsp = NULL, **SgmApsp = NULL, **SgmLpsp = NULL, **SgmGpsp = NULL;
matfunc_dcmplx **GrnRpsp = NULL, **GrnApsp = NULL, **GrnLpsp = NULL, **GrnGpsp = NULL;
matfunc_dcmplx **SgmRqus = NULL, **SgmAqus = NULL, **SgmLqus = NULL, **SgmGqus = NULL;
matfunc_dcmplx **GrnRqus = NULL, **GrnAqus = NULL, **GrnLqus = NULL, **GrnGqus = NULL;
matfunc_dcmplx **SgmRqusPT = NULL, **SgmAqusPT = NULL, **SgmLqusPT = NULL, **SgmGqusPT = NULL;
matfunc_dcmplx **GrnRqusPT = NULL, **GrnAqusPT = NULL, **GrnLqusPT = NULL, **GrnGqusPT = NULL;
matfunc_dcmplx **LamR = NULL, **LamA = NULL, **LamL = NULL, **LamG = NULL;

static const char *varlabel = "variables:      ";
static const char *hybpspin = "hybpsp.in";
static const char *hmatin   = "hmat.in";
static const char *hpspin   = "hpsp.in";
static const char *xipspin  = "xipsp.in";
static const char *opticin  = "optic.in";

static void variables_alloc_sgm(void)
{
  // Pseudoparticle Self energies
  SgmRpsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(SgmRpsp, "SgmRpsp");
  SgmApsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(SgmApsp, "SgmApsp");
  SgmLpsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(SgmLpsp, "SgmLpsp");
  SgmGpsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(SgmGpsp, "SgmGpsp");
  // Pseudoparticle transformed quasiparticle Self energies
  SgmRqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmRqus, "SgmRqus");
  SgmAqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmAqus, "SgmAqus");
  SgmLqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmLqus, "SgmLqus");
  SgmGqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmGqus, "SgmGqus");
  // Perturbation theory quasiparticle Self energies
  SgmRqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmRqusPT, "SgmRqusPT");
  SgmAqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmAqusPT, "SgmAqusPT");
  SgmLqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmLqusPT, "SgmLqusPT");
  SgmGqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(SgmGqusPT, "SgmGqusPT");
  // Hybridization
  LamR = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(LamR, "LamR");
  LamA = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(LamA, "LamA");
  LamL = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(LamL, "LamL");
  LamG = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(LamG, "LamG");
  // Light coupling
  if(isOptic) {
    PiR = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    check_mem(PiR, "PiR");
    PiA = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    check_mem(PiA, "PiA");
    PiL = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    check_mem(PiL, "PiL");
    PiG = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
    check_mem(PiG, "PiG");
  }
  
  return;
  
 error:
  abort();
}

static void variables_alloc_grn(void)
{
  // Pseudoparticle Green functions
  GrnRpsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(GrnRpsp, "GrnRpsp");
  GrnApsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(GrnApsp, "GrnApsp");
  GrnLpsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(GrnLpsp, "GrnLpsp");
  GrnGpsp = (matfunc_dcmplx **) calloc(Nblk, sizeof(matfunc_dcmplx *));
  check_mem(GrnGpsp, "GrnGpsp");
  // Pseudoparticle transformed quasiparticle Self energies
  GrnRqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnRqus, "GrnRqus");
  GrnAqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnAqus, "GrnAqus");
  GrnLqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnLqus, "GrnLqus");
  GrnGqus = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnGqus, "GrnGqus");
  // Perturbation theory quasiparticle Self energies
  GrnRqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnRqusPT, "GrnRqusPT");
  GrnAqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnAqusPT, "GrnAqusPT");
  GrnLqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnLqusPT, "GrnLqusPT");
  GrnGqusPT = (matfunc_dcmplx **) calloc(Nspin, sizeof(matfunc_dcmplx *));
  check_mem(GrnGqusPT, "GrnGqusPT");
  
  return;
  
 error:
  abort();
}

static void variables_deall_sgm(void)
{
  int ispin, iblk;
  
  if(isOptic) {
    if(PiG) {
      for(ispin = Nspin - 1; ispin > -1; --ispin)
        if(PiG[ispin]) matfunc_del_dcmplx(PiG[ispin]);
      freeup(PiG);
    }
    if(PiL) {
      for(ispin = Nspin - 1; ispin > -1; --ispin)
        if(PiL[ispin]) matfunc_del_dcmplx(PiL[ispin]);
      freeup(PiL);
    }
    if(PiA) {
      for(ispin = Nspin - 1; ispin > -1; --ispin)
        if(PiA[ispin]) matfunc_del_dcmplx(PiA[ispin]);
      freeup(PiA);
    }
    if(PiR) {
      for(ispin = Nspin - 1; ispin > -1; --ispin)
        if(PiR[ispin]) matfunc_del_dcmplx(PiR[ispin]);
      freeup(PiR);
    }
  }
  
  if(LamG) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamG[ispin]);
    }
    freeup(LamG);
  }
  if(LamL) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamL[ispin]);
    }
    freeup(LamL);
  }
  if(LamA) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamA[ispin]);
    }
    freeup(LamA);
  }
  if(LamR) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(LamR[ispin]);
    }
    freeup(LamR);
  }
  
  if(SgmGqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmGqusPT[ispin]);
    }
    freeup(SgmGqusPT);
  }
  if(SgmLqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmLqusPT[ispin]);
    }
    freeup(SgmLqusPT);
  }
  if(SgmAqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmAqusPT[ispin]);
    }
    freeup(SgmAqusPT);
  }
  if(SgmRqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmRqusPT[ispin]);
    }
    freeup(SgmRqusPT);
  }
  
  if(SgmGqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmGqus[ispin]);
    }
    freeup(SgmGqus);
  }
  if(SgmLqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmLqus[ispin]);
    }
    freeup(SgmLqus);
  }
  if(SgmAqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmAqus[ispin]);
    }
    freeup(SgmAqus);
  }
  if(SgmRqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(SgmRqus[ispin]);
    }
    freeup(SgmRqus);
  }
  
  if(SgmGpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(SgmGpsp[iblk]);
    }
    freeup(SgmGpsp);
  }
  if(SgmLpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(SgmLpsp[iblk]);
    }
    freeup(SgmLpsp);
  }
  if(SgmApsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(SgmApsp[iblk]);
    }
    freeup(SgmApsp);
  }
  if(SgmRpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(SgmRpsp[iblk]);
    }
    freeup(SgmRpsp);
  }
}

static void variables_deall_grn(void)
{
  int ispin, iblk;
  
  if(GrnGqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnGqusPT[ispin]);
    }
    freeup(GrnGqusPT);
  }
  if(GrnLqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnLqusPT[ispin]);
    }
    freeup(GrnLqusPT);
  }
  if(GrnAqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnAqusPT[ispin]);
    }
    freeup(GrnAqusPT);
  }
  if(GrnRqusPT) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnRqusPT[ispin]);
    }
    freeup(GrnRqusPT);
  }
  
  if(GrnGqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnGqus[ispin]);
    }
    freeup(GrnGqus);
  }
  if(GrnLqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnLqus[ispin]);
    }
    freeup(GrnLqus);
  }
  if(GrnAqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnAqus[ispin]);
    }
    freeup(GrnAqus);
  }
  if(GrnRqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      matfunc_del_dcmplx(GrnRqus[ispin]);
    }
    freeup(GrnRqus);
  }
  
  if(GrnGpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(GrnGpsp[iblk]);
    }
    freeup(GrnGpsp);
  }
  if(GrnLpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(GrnLpsp[iblk]);
    }
    freeup(GrnLpsp);
  }
  if(GrnApsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(GrnApsp[iblk]);
    }
    freeup(GrnApsp);
  }
  if(GrnRpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      matfunc_del_dcmplx(GrnRpsp[iblk]);
    }
    freeup(GrnRpsp);
  }
}

void variables_alloc(void)
{
  FILE         *fp = NULL;
  char         line[SHRT_MAX];
  int          ienv, ipsp, jpsp, ipspblk, jpspblk, iq, idim, jdim, iblk, jblk, ispin,
               blkpspi, npspblki, nn, ii;
  dreal        tmpre, tmpim;
  mat2d_dcmplx *hpspi = NULL, *xipspii = NULL, **xipspi = NULL, *chipspi = NULL;
  
  if(Rank == Root) {
    fp = fopen(hybpspin, "r");
    
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", &Efrm);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", &Temp);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf %lf %d", &Emin, &Emax, &Nmesh);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d",  &Nenv);
    
    Muenv = (dreal *) calloc(Nenv, sizeof(dreal));
    
    for(ienv = 0; ienv < Nenv; ++ienv) {
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", Muenv + ienv);
    }
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d %d", &Nspin, &Ndim);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d %d", &Nqus, &Nblk);
    
    Indqus = (int *) calloc(Nqus, sizeof(int));
    
    for(iq = 0; iq < Nqus; ++iq) {
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%d", Indqus + iq);
    }
    
    Npsp = numeric_ipow(2, Nspin*Nqus);
    Zeta = (int *) calloc(Npsp, sizeof(int));
    Blkpsp = (int *) calloc(Npsp, sizeof(int));
    
    for(ipsp = 0; ipsp < Npsp; ++ipsp) {
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%d %d", Zeta + ipsp, Blkpsp + ipsp);
    }
    
    Gamma = (dreal **) calloc(Nenv, sizeof(dreal *));
    for(ienv = 0; ienv < Nenv; ++ienv) {
      Gamma[ienv] = (dreal *) calloc(Ndim, sizeof(dreal));
      for(idim = 0; idim < Ndim; ++idim) {
        if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", Gamma[ienv] + idim);
      }
    }
    
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", &Mixing);
    
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d", &isOptic);
    
    fclose(fp);
  }
  
  sys_bcast_dreal(&Efrm,  1, Root);
  sys_bcast_dreal(&Temp,  1, Root);
  sys_bcast_dreal(&Emin,  1, Root);
  sys_bcast_dreal(&Emax,  1, Root);
  sys_bcast_dreal(&Mixing,1, Root);
  sys_bcast_int(&Nmesh,   1, Root);
  sys_bcast_int(&Nenv,    1, Root);
  sys_bcast_int(&Nspin,   1, Root);
  sys_bcast_int(&Ndim,    1, Root);
  sys_bcast_int(&Nqus,    1, Root);
  sys_bcast_int(&Nblk,    1, Root);
  sys_bcast_int(&Npsp,    1, Root);
  sys_bcast_int(&isOptic, 1, Root);
  
  check(Temp > 0.0, "Invalid Temp: %15.5e\n", Temp);
  check(Emax > Emin, "Invalid Emin/Emax: %15.5e/%15.5e\n", Emin, Emax);
  check(fabs(Emax) == fabs(Emin), "Invalid Emin/Emax: %15.5e/%15.5e\n", Emin, Emax);
  check(Nmesh > 0, "Invalid Nmesh: %5d\n", Nmesh);
  check(Nenv > 0, "Invalid Nenv: %5d\n", Nenv);
  check(Nspin > 0 && Nspin < 3, "Invalid Nspin: %5d\n", Nspin);
  check(Nqus > 0, "Invalid Nqus: %5d\n", Nqus);
  check(Nblk > 0 && Nblk <= Npsp, "Invalid Nblk: %5d", Nblk);
  
  if(Rank != Root) {
    Muenv = (dreal *) calloc(Nenv, sizeof(dreal));
    Zeta = (int *) calloc(Npsp, sizeof(int));
    Blkpsp = (int *) calloc(Npsp, sizeof(int));
    Indqus = (int *) calloc(Nqus, sizeof(int));
    Gamma = (dreal **) calloc(Nenv, sizeof(dreal *));
    for(ienv = 0; ienv < Nenv; ++ienv) {
      Gamma[ienv] = (dreal *) calloc(Ndim, sizeof(dreal));
    }
  }
  
  sys_bcast_dreal(Muenv, Nenv, Root);
  sys_bcast_int(Zeta, Npsp, Root);
  sys_bcast_int(Blkpsp, Npsp, Root);
  sys_bcast_int(Indqus, Nqus, Root);
  for(ienv = 0; ienv < Nenv; ++ienv) {
    sys_bcast_dreal(Gamma[ienv], Ndim, Root);
  }
  
  Efrm /= EV_HA;
  Temp *= HA_KBT;
  Beta  = 1.0 / Temp;
  Emin /= EV_HA;
  Emax /= EV_HA;
  for(ienv = 0; ienv < Nenv; ++ienv) {
    Muenv[ienv] /= EV_HA;
    for(idim = 0; idim < Ndim; ++idim) {
      Gamma[ienv][idim] /= EV_HA;
    }
  }
  
  Nmesh_aug = NAUG * (Nmesh - 1) + 1;
  Emax_aug  = NAUG * Emax;
  Emin_aug  = NAUG * Emin;
  
  Npspblk = (int *) calloc(Nblk, sizeof(int));
  check_mem(Npspblk, "Npspblk");
  for(ipsp = 0; ipsp < Npsp; ++ipsp) {
    blkpspi = Blkpsp[ipsp];
    ++Npspblk[blkpspi];
  }
  
  Pspblk = (int **) calloc(Nblk, sizeof(int *));
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki = Npspblk[iblk];
    Pspblk[iblk] = (int *) calloc(npspblki, sizeof(int));
  }
  for(iblk = 0; iblk < Nblk; ++iblk) {
    for(ipsp = 0, ipspblk = 0; ipsp < Npsp; ++ipsp) {
      jblk = Blkpsp[ipsp];
      if(jblk == iblk) {
        Pspblk[iblk][ipspblk] = ipsp;
        ++ipspblk;
      }
    }
  }
  
  // DMFT quasiparticle Hamiltonian
  Hmat = mat2d_new_dcmplx();
  check_mem(Hmat, "Hmat");
  Hmat->alloc(Ndim, Ndim, Hmat);
  
  if(Rank == Root) {
    fp = fopen(hmatin, "r");
    for(idim = 0; idim < Ndim; ++idim) {
      for(jdim = 0; jdim < Ndim; ++jdim) {
        if(fgets(line, sizeof(line), fp)) {
          sscanf(line, "%lf %lf", &tmpre, &tmpim);
        }
        Hmat->ptr[idim][jdim] = tmpre + I * tmpim;
      }
    }
    fclose(fp);
  }
  
  sys_bcast_dcmplx(Hmat->addr, Ndim*Ndim, Root);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(jdim = 0; jdim < Ndim; ++jdim) {
      Hmat->ptr[idim][jdim] /= EV_HA;
    }
  }
  
  // Pseudoparticle Hamiltonian
  Hpsp = (mat2d_dcmplx **) calloc(Nblk, sizeof(mat2d_dcmplx *));
  check_mem(Hpsp, "Hpsp");
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki = Npspblk[iblk];
    hpspi = mat2d_new_dcmplx();
    check_mem(hpspi, "hpspi");
    hpspi->alloc(npspblki, npspblki, hpspi);
    Hpsp[iblk] = hpspi;
  }
  
  if(Rank == Root) {
    fp = fopen(hpspin, "r");
    
    for(iblk = 0; iblk < Nblk; ++iblk) {
      npspblki = Npspblk[iblk];
      hpspi = Hpsp[iblk];
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
          if(fgets(line, sizeof(line), fp)) {
            sscanf(line, "%lf %lf", &tmpre, &tmpim);
          }
          hpspi->ptr[ipspblk][jpspblk] = tmpre + I * tmpim;
        }
      }
    }
    
    fclose(fp);
  }
  
  Epsp = (dreal *) calloc(Npsp, sizeof(dreal));
  
  for(iblk = 0; iblk < Nblk; ++iblk) {
    npspblki = Npspblk[iblk];
    hpspi = Hpsp[iblk];
    sys_bcast_dcmplx(hpspi->addr, npspblki*npspblki, Root);
    for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
      for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
        hpspi->ptr[ipspblk][jpspblk] /= EV_HA;
      }
      ipsp = Pspblk[iblk][ipspblk];
      Epsp[ipsp] = (dreal) hpspi->ptr[ipspblk][ipspblk];
    }
  }
  
  // Pseudoparticle overlap matrix
  Xipsp = (mat2d_dcmplx ***) calloc(Nspin, sizeof(mat2d_dcmplx **));
  check_mem(Xipsp, "Xipsp");
  for(ispin = 0; ispin < Nspin; ++ispin) {
    xipspi = (mat2d_dcmplx **) calloc(Nqus, sizeof(mat2d_dcmplx *));
    check_mem(xipspi, "xipspi");
    for(iq = 0; iq < Nqus; ++iq) {
      xipspii = mat2d_new_dcmplx();
      check_mem(xipspii, "xipspii");
      xipspii->alloc(Npsp, Npsp, xipspii);
      xipspi[iq] = xipspii;
    }
    Xipsp[ispin] = xipspi;
  }
  
  if(Rank == Root) {
    fp = fopen(xipspin, "r");
    
    for(ispin = 0; ispin < Nspin; ++ispin) {
      xipspi = Xipsp[ispin];
      for(iq = 0; iq < Nqus; ++iq) {
        xipspii = xipspi[iq];
        if(fgets(line, sizeof(line), fp)) {
          sscanf(line, "%d", &nn);
        }
        for(ii = 0; ii < nn; ++ii) {
          if(fgets(line, sizeof(line), fp)) {
            sscanf(line, "%d %d %lf %lf", &jpsp, &ipsp, &tmpre, &tmpim);
          }
          xipspii->ptr[ipsp][jpsp] = tmpre + I * tmpim;
        }
      }
    }
    
    fclose(fp);
  }
  
  for(ispin = 0; ispin < Nspin; ++ispin) {
    xipspi = Xipsp[ispin];
    for(iq = 0; iq < Nqus; ++iq) {
      xipspii = xipspi[iq];
      sys_bcast_dcmplx(xipspii->addr, Npsp*Npsp, Root);
    }
  }
  
  Sumdospsp = (dreal *) calloc(Npsp, sizeof(dreal));
  check_mem(Sumdospsp, "Sumdospsp");
  Dospsp = (dreal **) calloc(Npsp, sizeof(dreal *));
  check_mem(Dospsp, "Dospsp");
  
  Sumpoppsp = (dreal *) calloc(Npsp, sizeof(dreal));
  check_mem(Sumpoppsp, "Sumpoppsp");
  Poppsp = (dreal **) calloc(Npsp, sizeof(dreal *));
  check_mem(Poppsp, "Poppsp");
  
  Dosqus = (dreal ***) calloc(Nspin, sizeof(dreal **));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    Dosqus[ispin] = (dreal **) calloc(Nqus, sizeof(dreal *));
  }
  Popqus = (dreal ***) calloc(Nspin, sizeof(dreal **));
  for(ispin = 0; ispin < Nspin; ++ispin) {
    Popqus[ispin] = (dreal **) calloc(Nqus, sizeof(dreal *));
  }
  
  Itot = (dcmplx *) calloc(Nenv, sizeof(dcmplx));
  
  if(isOptic) {
    if(Rank == Root) {
      fp = fopen(opticin, "r");
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf", &TempOptic);
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%d %d", &Gstat, &Xstat);
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf %lf", &OmegOptic, &RateOptic);
      fclose(fp);
    }
    sys_bcast_int(&Gstat, 1, Root);
    sys_bcast_int(&Xstat, 1, Root);
    sys_bcast_dreal(&TempOptic, 1, Root);
    sys_bcast_dreal(&OmegOptic, 1, Root);
    sys_bcast_dreal(&RateOptic, 1, Root);
    
    TempOptic *= HA_KBT;
    BetaOptic  = 1.0 / TempOptic;
    OmegOptic /= EV_HA;
    RateOptic /= EV_HA;
    
    Chipsp = (mat2d_dcmplx **) calloc(Nspin, sizeof(mat2d_dcmplx *));
    for(ispin = 0; ispin < Nspin; ++ispin) {
      chipspi = mat2d_new_dcmplx();
      chipspi->alloc(Npsp, Npsp, chipspi);
      lapack_zgemm(Npsp, Npsp, Npsp, 'N', 'C', 1.0, Xipsp[ispin][Xstat]->addr, Xipsp[ispin][Gstat]->addr,
                   0.0, chipspi->addr);
      Chipsp[ispin] = chipspi;
    }
  }
  
  variables_alloc_sgm();
  variables_alloc_grn();
  
  return;
  
 error:
  abort();
}

void variables_deall(void)
{
  int ipsp, iblk, ienv, ispin, iq;
  
  variables_deall_grn();
  variables_deall_sgm();
  
  if(isOptic) {
    if(Chipsp) {
      for(ispin = Nspin - 1; ispin > -1; --ispin)
        mat2d_del_dcmplx(Chipsp[ispin]);
      freeup(Chipsp);
    }
  }
  
  if(Popqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      if(Popqus[ispin]) {
        for(iq = Nqus - 1; iq > -1; --iq) {
          if(Popqus[ispin][iq]) freeup(Popqus[ispin][iq]);
        }
        freeup(Popqus[ispin]);
      }
    }
    freeup(Popqus);
  }
  if(Dosqus) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      if(Dosqus[ispin]) {
        for(iq = Nqus - 1; iq > -1; --iq) {
          if(Dosqus[ispin][iq]) freeup(Dosqus[ispin][iq]);
        }
        freeup(Dosqus[ispin]);
      }
    }
    freeup(Dosqus);
  }
  if(Poppsp) {
    for(ipsp = Npsp - 1; ipsp > -1; --ipsp) {
      if(Poppsp[ipsp]) freeup(Poppsp[ipsp]);
    }
    freeup(Poppsp);
  }
  if(Dospsp) {
    for(ipsp = Npsp - 1; ipsp > -1; --ipsp) {
      if(Dospsp[ipsp]) freeup(Dospsp[ipsp]);
    }
    freeup(Dospsp);
  }
  if(Xipsp) {
    for(ispin = Nspin - 1; ispin > -1; --ispin) {
      for(iq = Nqus - 1; iq > -1; --iq) {
        if(Xipsp[ispin][iq]) mat2d_del_dcmplx(Xipsp[ispin][iq]);
      }
      freeup(Xipsp[ispin]);
    }
    freeup(Xipsp);
  }
  if(Hpsp) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      if(Hpsp[iblk]) mat2d_del_dcmplx(Hpsp[iblk]);
    }
    freeup(Hpsp);
  }
  if(Pspblk) {
    for(iblk = Nblk - 1; iblk > -1; --iblk) {
      if(Pspblk[iblk]) freeup(Pspblk[iblk]);
    }
    freeup(Pspblk);
  }
  if(Gamma) {
    for(ienv = Nenv - 1; ienv > -1; --ienv) {
      if(Gamma[ienv]) freeup(Gamma[ienv]);
    }
    freeup(Gamma);
  }
  
  mat2d_del_dcmplx(Hmat);
  
  if(Itot)      freeup(Itot);
  if(Sumpoppsp) freeup(Sumpoppsp);
  if(Sumdospsp) freeup(Sumdospsp);
  if(Epsp)      freeup(Epsp);
  if(Npspblk)   freeup(Npspblk);
  if(Blkpsp)    freeup(Blkpsp);
  if(Zeta)      freeup(Zeta);
  if(Indqus)    freeup(Indqus);
  if(Muenv)     freeup(Muenv);
}

void variables_print(void)
{
  FILE         *fp = NULL;
  int          ienv, ipsp, jpsp, ipspblk, jpspblk, iq, idim, iblk, ispin, npspblki;
  dreal        tmpre, tmpim;
  char         *hpspname   = "Hpsp.out";
  char         *xipspname  = "Xipsp.out";
  char         *chipspname = "Chipsp.out";
  char         *hmatname   = "Hmat.out";
  mat2d_dcmplx *hpspi = NULL, *xipspii = NULL, **xipspi = NULL, *chipspi = NULL;
  
  if(Rank == Root) {
    
    printf("%s\n", varlabel);
    printf("%s Efrm (Ha)    = %15.5e\n", varlabel, Efrm);
    printf("%s Temp (Ha)    = %15.5e\n", varlabel, Temp);
    printf("%s Beta         = %15.5e\n", varlabel, Beta);
    printf("%s Emin (Ha)    = %15.5e\n", varlabel, Emin);
    printf("%s Emax (Ha)    = %15.5e\n", varlabel, Emax);
    printf("%s Nmesh        = %5d\n",    varlabel, Nmesh);
    printf("%s Emin_aug (Ha)= %15.5e\n", varlabel, Emin_aug);
    printf("%s Emax_aug (Ha)= %15.5e\n", varlabel, Emax_aug);
    printf("%s Nmesh_aug    = %5d\n",    varlabel, Nmesh_aug);
    printf("%s Nenv         = %5d\n",    varlabel, Nenv);
    printf("%s Muenv (Ha)   =\n",        varlabel);
    for(ienv = 0; ienv < Nenv; ++ienv) {
      printf("%s %5d %15.5e\n", varlabel, ienv, Muenv[ienv]);
    }
    printf("%s Nspin        = %5d\n",    varlabel, Nspin);
    printf("%s Ndim         = %5d\n",    varlabel, Ndim);
    printf("%s Nqus         = %5d\n",    varlabel, Nqus);
    printf("%s Npsp         = %5d\n",    varlabel, Npsp);
    printf("%s Indqus       =\n",        varlabel);
    for(iq = 0; iq < Nqus; ++iq) {
      printf("%s %5d %5d\n", varlabel, iq, Indqus[iq]);
    }
    printf("%s Zeta, Blkpsp =\n",        varlabel);   
    for(ipsp = 0; ipsp < Npsp; ++ipsp) {
      printf("%s %5d %5d %5d\n", varlabel, ipsp, Zeta[ipsp], Blkpsp[ipsp]);
    }
    printf("%s Npspblk      =\n",        varlabel);
    for(iblk = 0; iblk < Nblk; ++iblk) {
      printf("%s %5d %5d\n", varlabel, iblk, Npspblk[iblk]);
    }
    printf("%s Pspblk       =\n",        varlabel);
    for(iblk = 0; iblk < Nblk; ++iblk) {
      npspblki = Npspblk[iblk];
      printf("%s %5d ", varlabel, iblk);
      for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
        printf("%5d ", Pspblk[iblk][ipspblk]);
      }
      printf("\n");
    }
    printf("%s Epsp (Ha)    =\n",        varlabel);
    for(ipsp = 0; ipsp < Npsp; ++ipsp) {
      printf("%s %5d %15.5e\n", varlabel, ipsp, Epsp[ipsp]);
    }
    printf("%s Gamma (Ha)   =\n",        varlabel);
    for(ienv = 0; ienv < Nenv; ++ienv) {
      for(idim = 0; idim < Ndim; ++idim) {
        printf("%s %5d %5d %15.5e\n", varlabel, ienv, idim, Gamma[ienv][idim]);
      }
    }
    printf("%s Mixing       = %15.5e\n", varlabel, Mixing);
    printf("%s\n", varlabel);
    if(isOptic) {
      printf("%s TempOptic (Ha) = %15.5e\n",  varlabel, TempOptic);
      printf("%s BetaOptic (Ha) = %15.5e\n",  varlabel, BetaOptic);
      printf("%s OmegOptic (Ha) = %15.5e\n",  varlabel, OmegOptic);
      printf("%s RateOptic (Ha) = %15.5e\n",  varlabel, RateOptic);
      printf("%s Gstat, Xstat   = %5d %5d\n", varlabel, Gstat, Xstat);
      printf("%s\n", varlabel);
    }
    fflush(stdout);
    // For debug
    fp = fopen(hpspname, "w");
    for(iblk = 0; iblk < Nblk; ++iblk) {
      hpspi = Hpsp[iblk];
      npspblki = Npspblk[iblk];
      fprintf(fp, "%s: %5d\n", "Block", iblk);
      for(jpspblk = 0; jpspblk < npspblki; ++jpspblk) {
        for(ipspblk = 0; ipspblk < npspblki; ++ipspblk) {
          tmpre = creal(hpspi->ptr[ipspblk][jpspblk]);
          tmpim = cimag(hpspi->ptr[ipspblk][jpspblk]);
          fprintf(fp, " (%15.5e %15.5e) ", tmpre, tmpim);
        }
        fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    
    fp = fopen(xipspname, "w");
    for(ispin = 0; ispin < Nspin; ++ispin) {
      xipspi = Xipsp[ispin];
      for(iq = 0; iq < Nqus; ++iq) {
        xipspii = xipspi[iq];     
        fprintf(fp, "%s %5d %5d\n", "Spin, Quas:", ispin, iq);
        for(jpsp = 0; jpsp < Npsp; ++jpsp) {
          for(ipsp = 0; ipsp < Npsp; ++ipsp) {
            tmpre = creal(xipspii->ptr[ipsp][jpsp]);
            tmpim = cimag(xipspii->ptr[ipsp][jpsp]);
            fprintf(fp, " (%15.5e %15.5e) ", tmpre, tmpim);
          }
          fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
      }
    }
    fclose(fp);
    
    if(isOptic) {
      fp = fopen(chipspname, "w");
      for(ispin = 0; ispin < Nspin; ++ispin) {
        chipspi = Chipsp[ispin];
        fprintf(fp, "%s %5d\n", "Spin:", ispin);
        for(jpsp = 0; jpsp < Npsp; ++jpsp) {
          for(ipsp = 0; ipsp < Npsp; ++ipsp) {
            tmpre = creal(chipspi->ptr[ipsp][jpsp]);
            tmpim = cimag(chipspi->ptr[ipsp][jpsp]);
            fprintf(fp, " (%15.5e %15.5e) ", tmpre, tmpim);
          }
          fprintf(fp, "\n");
        }
      }
      fclose(fp);
    }
  }
  
  Hmat->print(hmatname, Hmat);
}
