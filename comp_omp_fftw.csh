gcc -o $1_omp $1.c -std=c99 -fopenmp -D_OMPTHREAD_ -Wall -O3 -march=native -lfftw3_threads -lm
