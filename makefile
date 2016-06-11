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

aux = auxiliary.o Input_variables.o 


simfast21: get_densityfield get_velocityfield \
     get_halos get_nldensity adjust_halos get_halo_deltan\
     get_HIIbubbles get_SFR \
     xalpha xc epsilonXon integratexe \
     integrateTempX t21 opdepth get_dndm get_dndm_nbins adndm abias \
     power3d rz

%.o: %.c auxiliary.h Input_variables.h
	${cc} ${flags} -c -o $@ $<

get_densityfield: $(aux) get_densityfield.o
	$(cc) -o get_densityfield.x get_densityfield.o $(aux) $(flags)

get_velocityfield: $(aux) get_velocityfield.o
	$(cc) -o get_velocityfield.x get_velocityfield.o $(aux) $(flags)

get_halos: $(aux) get_halos.o
	$(cc) -o get_halos.x get_halos.o $(aux) $(flags)

adjust_halos: $(aux) adjust_halos.o
	$(cc) -o adjust_halos.x adjust_halos.o $(aux) $(flags)

get_nldensity: $(aux) get_nldensity.o
	$(cc) -o get_nldensity.x get_nldensity.o $(aux) $(flags)

get_HIIbubbles: $(aux) get_HIIbubbles.o user_functions.o
	$(cc) -o get_HIIbubbles.x get_HIIbubbles.o user_functions.o $(aux) $(flags)

get_SFR: $(aux) get_SFR.o
	$(cc) -o get_SFR.x get_SFR.o $(aux) $(flags)

xalpha: $(aux) xalpha.o
	$(cc) -o xalpha.x xalpha.o $(aux) $(flags)

xc: $(aux) xc.o
	$(cc) -o xc.x xc.o $(aux) $(flags)

epsilonXon: $(aux) epsilonXon.o
	$(cc) -o epsilonXon.x epsilonXon.o $(aux) $(flags)

integratexe: $(aux) integratexe.o
	$(cc) -o integratexe.x integratexe.o $(aux) $(flags)

integrateTempX: $(aux) integrateTempX.o
	$(cc) -o integrateTempX.x integrateTempX.o $(aux) $(flags)

t21: $(aux) t21.o
	$(cc) -o t21.x t21.o $(aux) $(flags)

power3d: tools/power3d.o
	$(cc) -o tools/power3d.x tools/power3d.o $(flags)


#tools
get_halo_deltan: $(aux) tools/get_halo_deltan.o
	$(cc) -o tools/get_halo_deltan.x tools/get_halo_deltan.o $(aux) $(flags)

get_dndm: $(aux) tools/get_dndm.o
	$(cc) -o tools/get_dndm.x tools/get_dndm.o $(aux) $(flags)

get_dndm_nbins: $(aux) tools/get_dndm_nbins.o
	$(cc) -o tools/get_dndm_nbins.x tools/get_dndm_nbins.o $(aux) $(flags)

adndm: $(aux) tools/adndm.o
	$(cc) -o tools/adndm.x tools/adndm.o $(aux) $(flags)

abias: $(aux) tools/abias.o
	$(cc) -o tools/abias.x tools/abias.o $(aux) $(flags)

rz: $(aux) tools/rz.o
	$(cc) -o tools/rz.x tools/rz.o $(aux) $(flags)

power3d: tools/power3d.o
	$(cc) -o tools/power3d.x tools/power3d.o $(flags)

opdepth: Input_variables.o tools/opdepth.o
	$(cc) -o tools/opdepth.x tools/opdepth.o Input_variables.o -std=c99 -Wall -O3 -lm 

clean:
	rm *.o *.x tools/*.o tools/*.x


