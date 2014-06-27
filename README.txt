
*************************************************************************************************************************
*                                                                                                                       *
*  SimFast21 v1.0                                                                                                       *
*                                                                                                                       *   
*  Copyright 2011 Mario Santos, Luis Ferramacho, Marta Silva, Alexandre Amblard, Asantha Cooray                         *
*                                                                                                                       *
*  This program is distributed under the terms of the GNU General Public License (http://www.gnu.org/copyleft/gpl.html) *
*                                                                                                                       *
*  Contact: mgrsantos@ist.utl.pt     www.simfast21.org                                                                  *
*                                                                                                                       *
*************************************************************************************************************************


= Description =

This program generates a simulation of the cosmological 21cm signal (see http://arxiv.org/abs/0911.2219 for details).
To run just do: ./simfast21 base_dir where "base_dir" is the base directory where you want your simulation to reside (note that the files can be quite big). This base directory must contain the parameters file simfast21.ini (see the copy that comes with the distribution for further details on how to use this file).
The program will create the following directory structure:
base_dir/density     	    - contains the matter density (delta) box for z=0 and nonlinear boxes for all redshifts
        /Velocity    	    - Velocity at z=0
        /Halos       	    - halo catalogs and nonlinear collapsed mass boxes
        /Ionization  	    - ionization fraction boxes
	/SFR         	    - Star formation rate density boxes
	/xrays       	    - boxes to calculate heating due to xrays
	/Lya	     	    - Lya coupling
	/xc		    - collisional coupling
	/deltaTb       	    - final 21cm brightness temperature boxes
	/Output_text_files  - Several text files generated during the simulation (averages, etc)



= Compilation =

Just run make.
It requires the Gnu Scientific Library and the FFTW libraries.


