#! /bin/bash
mpirun -n $1 --bynode -hostfile hosts ./bin/simwave -v --grid 100,100,100 --iter 10  --epsilon 1e-17
