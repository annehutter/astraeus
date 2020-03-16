# astraeus
Home of the semi-analytical / semi-numerical galaxy evolution and reionization code ASTRAEUS (semi-numerical r**A**diative tran**S**fer coupling of galaxy forma**T**ion and **R**eionization in N-body d**A**rk matt**E**r sim**U**lation**S**)

When you should use this code
=============================

If you want to compute the self-consistent evolution of high-redshift galaxy properties and reionization, i.e. the galaxy properties in the presence of an inhomogeneous ultraviolet background and when and how (HI, HeII, HeIII) reionization occurs, then you should use this code. You will need 

- merger trees generated with ``cuteNresort``
- (optional) cosmological box with DM/gas overdensities **or** gas densities (grid); if none are given, DM/gas density is assumed to be homogeneous.

Why should you use it
=====================

1. **Modular** The code is written modular fashion, i.e. different modules can be switched on or off or chosen (see below for the inifiles).
2. **MPI Parallel** The code can be run on multiple cores and distributed memory.

Installation
============

Pre-requisities
---------------

Serial run
``````````

1. fftw3 library: ``fftw3 >= 3.3.3``
2. gsl library (math, integrate): ``gsl >= 1.16``

Parallel run
````````````

1. MPI library
2. fftw3 & fftw3-mpi library: ``fftw3 >= 3.3.3``
3. gsl library (math, integrate): ``gsl >= 1.16``

FFTW3
'''''

Go to the `FFTW webpage <http://www.fftw.org/download.html>`__ to install fftw3. Ensure to compile the library with the ``enable-mpi`` flag for parallel runs
::
    
    $ ./configure --enable-mpi
    $ make
    $ make install
    
Note: To create the dynamic libraries, run configure with the ``--enable-shared`` flag. 
    
GSL
'''

Go to the `GSL webpage <https://www.gnu.org/software/gsl/>`__ to install gsl and follow the instructions there. 


Download & Build
----------------

::

    $ git clone https://github.com/annehutter/grid-model.git
    $ make

This will download the code and first test case from the github directory and compile the source code.

Execution
---------

The first test case can then be run by
::

    $ ./astraeus iniFile.ini

``iniFile.ini`` contains all input parameters that are needed for any runs. For a different simulation the code does not need to be recompiled but only the parameter file iniFile.ini to be adapted.
