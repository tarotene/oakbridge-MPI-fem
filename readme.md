Run the commands below on the root directory (:= where this markdown is located).

```sh:commands
./sh/compile.sh

# parallelization test
./sh/check_parallelized.sh # Check termination by pjstat after a moment. Run with a small size input for safe terminzation.

# parallelization check
tail -n +2 ./out/out-serial.dat | cut -d "," -f 2 > _phi-serial.dat
tail -n +2 ./out/out-parallel.dat | cut -d "," -f 3 > _phi-parallel.dat
diff _phi-serial.dat _phi-parallel.dat > ./out/out-diff.dat
rm _phi-*.dat

# performance measurement
./sh/check_perf.sh # Run with a big size input for good performance.
```

