#compiler
#cc = /usr/local/bin/gcc-8
cc = gcc

#complilation flags

# if FFTW uses openmp:
fftwf = -lfftw3f_omp -lfftw3f
fftwd = -lfftw3_omp -lfftw3
# if FFTW used pthreads:
#fftwf = -lfftw3f_threads -lfftw3f
#fftwd = -lfftw3_threads -lfftw3

omp = -fopenmp -D_OMPTHREAD_
#omp = -Xpreprocessor -fopenmp -D_OMPTHREAD_
gsl = -lgsl -lgslcblas -lm
deb = -g

#flags = -std=c99 -Wall -O3 -march=native $(fftwd) $(fftwf) $(omp) $(gsl) -lm
flags = -std=c99 -Wall -O3 $(omp) $(fftwd) $(fftwf) $(gsl) -I .

auxo = auxiliary.o Input_variables.o 
auxh = auxiliary.h Input_variables.h 

main = get_densityfield.x get_velocityfield.x get_halos.x get_nldensity.x adjust_halos.x \
	xalpha.x xc.x epsilonXon.x integratexe.x integrateTempX.x t21.x get_HIIbubbles.x get_SFR.x

tools = tools/get_halo_deltan.x tools/opdepth.x tools/power3d_f.x tools/power3d_cross.x \
	tools/get_dndm_nbins.x tools/adndm.x tools/abias.x tools/power3d.x tools/rz.x tools/vel_grad.x

Simfast21: $(main) $(tools)

get_HIIbubbles.x: get_HIIbubbles.o user_functions.o $(auxo)
	$(cc) -o get_HIIbubbles.x get_HIIbubbles.o user_functions.o $(auxo) $(flags)

get_SFR.x: get_SFR.o user_functions.o $(auxo)
	$(cc) -o get_SFR.x get_SFR.o user_functions.o $(auxo) $(flags)

xalpha.x: xalpha.o user_functions.o $(auxo)
	$(cc) -o xalpha.x xalpha.o user_functions.o $(auxo) $(flags)

%.x: %.o $(auxo)
	$(cc) -o $@ $< $(auxo) $(flags)

%.o: %.c $(auxh)
	${cc} ${flags} -c $< -o $@

.PRECIOUS: %.o

clean:
	rm *.o *.x tools/*.o tools/*.x


