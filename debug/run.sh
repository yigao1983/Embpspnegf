#!/bin/bash

#valgrind --leak-check=full $HOME/Works/hybpspnegf-optic/obj/hybpspnegfexec >& Log &
#$HOME/Works/hybpspnegf-optic/obj/hybpspnegfexec >& Log &
time mpiexec -np 8 $HOME/Works/hybpspnegf-optic/obj/hybpspnegfexec >& Log &

# end of run.sh
