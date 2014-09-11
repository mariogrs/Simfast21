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
     get_halos get_nldensity adjust_halos \
     get_HIIbubbles get_SFR \
     xalpha xc epsilonXon integratexe \
     integrateTempX t21 opdepth

%.o: %.c auxiliary.h Input_variables.h
	${cc} ${flags} -c -o $@ $<

get_densityfield: $(aux) get_densityfield.o
	$(cc) -o get_densityfield.x get_densityfield.o $(aux) $(flags)

get_velocityfield: $(aux) get_velocityfield.o
	$(cc) -o get_velocityfield.x get_velocityfield.o $(aux) $(flags)

get_halos: $(aux) get_halos_fcoll.o
	$(cc) -o get_halos.x get_halos_fcoll.o $(aux) $(flags)

adjust_halos: $(aux) adjust_halos.o
	$(cc) -o adjust_halos.x adjust_halos.o $(aux) $(flags)

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

opdepth: Input_variables.o opdepth.o
	$(cc) -o opdepth.x opdepth.o Input_variables.o -std=c99 -Wall -O3 -lm 

clean:
	rm *.o get_densityfield.x get_velocityfield.x \
                   get_halos.x \
                   get_nldensity.x adjust_halos.x \
                   get_HIIbubbles.x\
                   get_SFR.x \
                   xalpha.x xc.x epsilonXon.x integratexe.x \
                   integrateTempX.x t21.x opdepth.x

