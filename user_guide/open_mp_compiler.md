
A common set of environment variables for parallel with open MP is.

```bash
export OMP_NUM_THREADS=8 # or number of threads you have.
export FFLAGS="-O2 -fopenmp -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8 -fallow-argument-mismatch"

export SLURM_CPU_BIND_VERBOSE=True
export OMP_PROCESS_BIND=true
export OMP_PLACES=threads
```
