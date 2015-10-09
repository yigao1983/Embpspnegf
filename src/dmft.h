#ifndef dmft_h
#define dmft_h

#include "precision.h"

extern matfunc_dcmplx **PiR_aug, **PiA_aug, **PiL_aug, **PiG_aug;
extern matfunc_dcmplx **LamR_aug, **LamA_aug, **LamL_aug, **LamG_aug;

void dmft_matfunc_aug(const matfunc_dcmplx *, matfunc_dcmplx *);
void dmft_sumdosqus(matfunc_dcmplx * const *, dreal **);
void dmft_sumpopqus(matfunc_dcmplx * const *, dreal **);
void dmft_alloc(void);
void dmft_deall(void);
void dmft_main(void);

#endif
