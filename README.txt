
*************************************************************************************************************************
*                                                                                                                       *
*  SimFast21 - new version to match code used in arXiv:1510.04280v1                                                     *                                                
*  Initial versions based on algorithms described in: http://arxiv.org/abs/0911.2219 and http://arxiv.org/abs/0708.2424 *
*														     *
*  This program is distributed under the terms of the GNU General Public License (http://www.gnu.org/copyleft/gpl.html) *
*                                                                                                                       *
*  Contact: mariogrs@gmail.com    https://github.com/mariogrs/Simfast21                                                 *                 *														     *
*                                                                                                                       *
*************************************************************************************************************************

= Code Contributors =
Mario Santos
Luis Ferramacho
Marta Silva
Alexandre Amblard
Sultan Hassan


= Description =

This program generates a simulation of the cosmological 21cm signal during the epoch of reionization.
To run just do: "./simfast21 base_dir" where "base_dir" is the base directory where you want your simulation to reside (note that the files can be quite big). This base directory must contain the parameters file simfast21.ini (see the copy that comes with the distribution for further details on how to use this file).
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


The folder "tools/" contains extra tools that can be used to analyse
the output boxes.

= Compilation =

Just run make. You might need to edit "makefile" in order to compile in your system.
It requires the Gnu Scientific Library and the FFTW libraries. Make
sure you compile the FFTW with --enable-threads or --enable-openmp and
use the appropriate flags in the make file.

See "output_files.txt" and "code_algorithm_diagram.png" for more
details on the code.



