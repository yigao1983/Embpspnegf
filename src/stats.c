#include <stdlib.h>
#include <math.h>
#include "stats.h"

dreal stats_ferm(dreal beta, dreal mu, dreal eng)
{
  dreal deng, ferm;
  
  deng = eng - mu;
  
  ferm = 1.0 / (exp(beta*deng) + 1.0);
  
  return ferm;
}

dreal stats_invsferm(dreal beta, dreal mu, dreal eng)
{
  dreal ferm;
  
  ferm = stats_ferm(beta, mu, eng);
  
  return (1.0 - ferm);
}

dreal stats_bose(dreal beta, dreal mu, dreal eng)
{
  dreal deng, bose;
  
  deng = eng - mu;
  
  bose = 1.0 / (exp(beta*deng) - 1.0);
  
  return bose;
}

dreal stats_invsbose(dreal beta, dreal mu, dreal eng)
{
  dreal bose;
  
  bose = stats_bose(beta, mu, eng);
  
  return (1.0 + bose);
}
