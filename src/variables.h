#ifndef variables_h
#define variables_h

#include "precision.h"
#include "matrix.h"
#include "matfunc.h"

extern int   Nmesh, Nenv, Nspin, Ndim, Nqus, Npsp, Nblk, Nmesh_aug;
extern int   isOptic, Gstat, Xstat;
extern dreal Efrm, Temp, Beta, Emin, Emax, Emin_aug, Emax_aug, Mixing;
extern dreal TempOptic, BetaOptic, OmegOptic, RateOptic;

extern int            *Npspblk, *Zeta, *Blkpsp, *Indqus;
extern int            **Pspblk;
extern dreal          *Muenv, *Ebuf, *Epsp, *Sumdospsp, *Sumpoppsp;
extern dreal          **Gamma, **Dospsp, **Poppsp;
extern dreal          ***Dosqus, ***Popqus;
extern dcmplx         *Itot;
extern mat2d_dcmplx   *Hmat;
extern mat2d_dcmplx   **Hpsp;
extern mat2d_dcmplx   ***Xipsp;
extern mat2d_dcmplx   **Chipsp;
extern matfunc_dcmplx **PiR, **PiA, **PiL, **PiG;
extern matfunc_dcmplx **SgmRpsp, **SgmApsp, **SgmLpsp, **SgmGpsp;
extern matfunc_dcmplx **GrnRpsp, **GrnApsp, **GrnLpsp, **GrnGpsp;
extern matfunc_dcmplx **SgmRqus, **SgmAqus, **SgmLqus, **SgmGqus;
extern matfunc_dcmplx **GrnRqus, **GrnAqus, **GrnLqus, **GrnGqus;
extern matfunc_dcmplx **SgmRqusPT, **SgmAqusPT, **SgmLqusPT, **SgmGqusPT;
extern matfunc_dcmplx **GrnRqusPT, **GrnAqusPT, **GrnLqusPT, **GrnGqusPT;
extern matfunc_dcmplx **LamR, **LamA, **LamL, **LamG;

void variables_alloc(void);
void variables_deall(void);
void variables_print(void);

#endif
