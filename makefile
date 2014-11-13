#compiler
cc = gcc

#complilation flags
fftwf = -lfftw3f -lfftw3f_threads
fftwd = -lfftw3 -lfftw3_threads
omp = -fopenmp -D_OMPTHREAD_
gsl = -lgsl -lgslcblas
deb = -g
#flags = -std=c99 -Wall -O3 -march=native $(fftwd) $(fftwf) $(omp) $(gsl) -lm
flags = -std=c99 -Wall -O3 $(fftwd) $(fftwf) $(omp) $(gsl) -lm

aux = auxiliary.o Input_variables.o 


simfast21: get_densityfield get_velocityfield \
     get_halos get_nldensity adjust_halos get_halo_deltan\
     get_HIIbubbles get_SFR \
     xalpha xc epsilonXon integratexe \
     integrateTempX t21 opdepth get_dndm get_dndmb get_dndmc adndm abias \
     power3df power3dd convert_halos rz

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

get_halo_deltan: $(aux) get_halo_deltan.o
	$(cc) -o get_halo_deltan.x get_halo_deltan.o $(aux) $(flags)

convert_halos: $(aux) convert_halos.o
	$(cc) -o convert_halos.x convert_halos.o $(aux) $(flags)

get_dndm: $(aux) get_dndm.o
	$(cc) -o get_dndm.x get_dndm.o $(aux) $(flags)

get_dndmb: $(aux) get_dndmb.o
	$(cc) -o get_dndmb.x get_dndmb.o $(aux) $(flags)

get_dndmc: $(aux) get_dndmc.o
	$(cc) -o get_dndmc.x get_dndmc.o $(aux) $(flags)

adndm: $(aux) adndm.o
	$(cc) -o adndm.x adndm.o $(aux) $(flags)

abias: $(aux) abias.o
	$(cc) -o abias.x abias.o $(aux) $(flags)

rz: $(aux) rz.o
	$(cc) -o rz.x rz.o $(aux) $(flags)

get_nldensity: $(aux) get_nldensity.o
	$(cc) -o get_nldensity.x get_nldensity.o $(aux) $(flags)

get_HIIbubbles: $(aux) get_HIIbubbles.o
	$(cc) -o get_HIIbubbles.x get_HIIbubbles.o $(aux) $(flags)

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

power3df: power3df.o
	$(cc) -o power3df.x power3df.o $(flags)

power3dd: power3dd.o
	$(cc) -o power3dd.x power3dd.o $(flags)

opdepth: Input_variables.o opdepth.o
	$(cc) -o opdepth.x opdepth.o Input_variables.o -std=c99 -Wall -O3 -lm 

clean:
	rm *.o *.x

