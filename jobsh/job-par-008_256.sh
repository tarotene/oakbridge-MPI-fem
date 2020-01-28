#!/bin/sh
#PJM -N "test"
#PJM -L rscgrp=lecture
#PJM -L node=8
#PJM --mpi proc=256
#PJM -L elapse=00:05:00
#PJM -g gt37
#PJM -j
#PJM -e err
#PJM -o ./stat/stat-par-008_256.lst
mpiexec.hydra -n ${PJM_MPI_PROC} ./bin/bin-parallel-noout-a.out
mpiexec.hydra -n ${PJM_MPI_PROC} ./bin/bin-parallel-noout-b.out
export I_MPI_PIN_PROCESSOR_LIST=0-15,28-43
mpiexec.hydra -n ${PJM_MPI_PROC} ./bin/bin-parallel-noout-a.out
mpiexec.hydra -n ${PJM_MPI_PROC} ./bin/bin-parallel-noout-b.out