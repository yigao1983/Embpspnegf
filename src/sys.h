#ifndef sys_h
#define sys_h

#include "precision.h"

extern int Size, Rank, Root;

void sys_init(int *, char ***);
void sys_quit(void);
void sys_sync(void);
void sys_abort(void);
void sys_bcast_char(char *, int, int);
void sys_bcast_int(int *, int, int);
void sys_bcast_dreal(dreal *, int, int);
void sys_bcast_dcmplx(dcmplx *, int, int);
void sys_allreduce_sum_int(int *, int);
void sys_allreduce_sum_dreal(dreal *, int);
void sys_allreduce_sum_dcmplx(dcmplx *, int);

#define abort() { printf("[ABORT] (%s:%s) Rank = %d\n", __FILE__, __FUNCTION__, Rank); sys_abort(); }

#endif
