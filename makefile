#compiler
cc = gcc

#complilation flags
fftwf = -lfftw3f -lfftw3f_threads
fftwd = -lfftw3 -lfftw3_threads
omp = -fopenmp -D_OMPTHREAD_
gsl = -lgsl -lgslcblas
deb = -g

#flags = -std=c99 -Wall -O3 -march=native $(fftwd) $(fftwf) $(omp) $(gsl) -lm
flags = -std=c99 -Wall -O3 $(fftwd) $(fftwf) $(omp) $(gsl) -lm -I .
flags2 = -std=c99 -Wall -O3 -lm -I .

auxo = auxiliary.o Input_variables.o 
auxh = auxiliary.h Input_variables.h 

main = get_densityfield.x get_velocityfield.x get_halos.x get_nldensity.x adjust_halos.x \
	xalpha.x xc.x epsilonXon.x integratexe.x integrateTempX.x t21.x get_HIIbubbles.x get_SFR.x

tools = tools/get_halo_deltan.x tools/opdepth.x tools/get_dndm.x \
	tools/get_dndm_nbins.x tools/adndm.x tools/abias.x tools/power3d.x tools/rz.x

Simfast21: $(main) $(tools)

get_HIIbubbles.x: get_HIIbubbles.o user_functions.o $(auxo)
	$(cc) -o get_HIIbubbles.x get_HIIbubbles.o user_functions.o $(auxo) $(flags)

get_SFR.x: get_SFR.o user_functions.o $(auxo)
	$(cc) -o get_SFR.x get_SFR.o user_functions.o $(auxo) $(flags)

%.x: %.o $(auxo)
	$(cc) -o $@ $< $(auxo) $(flags)

%.o: %.c $(auxh)
	${cc} ${flags} -c $< -o $@

clean:
	rm *.o *.x tools/*.o tools/*.x


