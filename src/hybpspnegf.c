#include <stdlib.h>
#include "sys.h"
#include "variables.h"
#include "dmft.h"
#include "perturb.h"

int main(int argc, char *argv[])
{
  // System initialization
  sys_init(&argc, &argv);
  // Variables allocation
  variables_alloc();
  variables_print();
  perturb_alloc();
  perturb_print();
  // DMFT
  dmft_alloc();
  dmft_main();
  dmft_deall();
  // Variables deallocation
  perturb_deall();
  variables_deall();
  // System quit
  sys_quit();
  
  return EXIT_SUCCESS;
}
