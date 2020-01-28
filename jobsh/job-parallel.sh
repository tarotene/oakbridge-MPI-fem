#!/bin/sh
#PJM -N "test"
#PJM -L rscgrp=lecture
#PJM -L node=2
#PJM --mpi proc=56
#PJM -L elapse=00:05:00
#PJM -g gt37
#PJM -j
#PJM -e err
#PJM -o ./stat/stat-parallel.lst
export I_MPI_PIN_PROCESSOR_LIST=0-27
mpiexec.hydra -n ${PJM_MPI_PROC} ./bin/bin-parallel-b.out