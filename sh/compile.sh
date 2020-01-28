#!/bin/bash

icc -O3 ./src/src-serial.c -o ./bin/bin-serial.out
mpiicc -align -O3 -axCORE-AVX512 ./src/src-parallel.c -o ./bin/bin-parallel-a.out
mpiicc -O3 ./src/src-parallel.c -o ./bin/bin-parallel-b.out
mpiicc -align -O3 -axCORE-AVX512 ./src/src-parallel-noout.c -o ./bin/bin-parallel-noout-a.out
mpiicc -O3 ./src/src-parallel.c -o ./bin/bin-parallel-noout-b.out