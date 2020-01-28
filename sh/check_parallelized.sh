#!/bin/bash

pjsub ./jobsh/job-parallel.sh
./bin/bin-serial.out > ./stat/stat-serial.lst