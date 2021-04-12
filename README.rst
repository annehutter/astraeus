ASTRAEUS
========

Home of the semi-analytical / semi-numerical galaxy evolution and reionization code ASTRAEUS (semi-numerical rAdiative tranSfer coupling of galaxy formaTion and Reionization in n-body dArk mattEr simUlationS). The underlying model is described in `Hutter et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.tmp..610H/abstract>`__. Suggestions for improvements are very welcome and can be sent to `this email address <a.k.hutter@rug.nl>`_.

When you should use this code
=============================

If you want to compute the self-consistent evolution of high-redshift galaxy properties and reionization, i.e. the galaxy properties in the presence of an inhomogeneous ultraviolet background and the time evolution of the ionization (HI, HeII, HeIII) fields, then you should use this code. You will need 

- merger trees generated with ``Consistent-Trees`` and resorted to be locally horizontal with ``cutNresort``.
- (optional) cosmological box with DM/gas overdensities **or** gas densities (grid); if none are given, DM/gas density is assumed to be homogeneous.

Note that ASTRAEUS can also be run on analytical merger trees when reionization is switched off (see below).

Why should you use it
=====================

1. **Modular** The code is written modular fashion, i.e. different modules can be switched on or off or chosen (see below for the inifiles).
2. **MPI Parallel** The code can be run on multiple cores and distributed memory.

Source files
============

The astraeus directory constains the Makefile as well as the executable ``astraeus`` after the source code in the `src <https://github.com/annehutter/astraeus/src>`__ directory has been compiled.

The directory `analysis <https://github.com/annehutter/astraeus/analysis>`__ constains an analysis code that runs on the ASTRAEUS tree output files. 

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

    $ git clone https://github.com/annehutter/astraeus.git
    $ make

This will download the code and a test case from the github directory and compile the source code.

Common.mk
---------

Adjust the paths of the FFTW libraries if they are not enabled by default or uncomment those lines and remove FFTW flags and links. 

Note that the ``OUTPUTSNAPNUMBER`` needs to be adjusted to the number of snapshots of the merger trees, i.e. if the merger trees start at snap ``0`` and end at snap ``i``, then set ``OUTPUTSNAPNUMBER = i`` while ``endSnapshot = i+1`` in the parameter file (see below).

Execution
---------

The first test case can then be run by
::

    $ ./astraeus iniFile.ini

``iniFile.ini`` contains the input parameters that are needed for any run. For a different simulation the code does not need to be recompiled but only the parameter file ``iniFile.ini`` to be adapted. ``iniFile.ini`` contains the input parameters for the galaxy evolution model and cosmology; furthermore it needs to contain the path to the ``cifogIniFile.ini``. In ``cifogIniFile.ini`` the input parameters for `CIFOG <https://ui.adsabs.harvard.edu/abs/2018ascl.soft03002H/abstract>`__ (described in `Hutter (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.1549H/abstract>`__) are set, for more details and what each parameter means please visit the `CIFOG github page <https://github.com/annehutter/grid-model/>`__.

Parameter File
==============

**[Input]**

- ``numFiles`` = *[integer]* number of input tree files
- ``fileName`` = path to input tree files; note the code assumes all the trees have the same name and differ just by their suffix ``_i``

**[Simulation]**

- ``redshiftFile`` = path to the redshift file which contains 3 columns: snapshot number, redshift z, scale factor a
- ``startSnapshot`` = *[integer]* snapshot number where to start computing the evolution of galaxy properties (and reionization) [default: 0]
- ``endSnapshot`` = *[integer]* snapshot where to end computing the evolution of galaxy properties (and reionization) [default: 74]
- ``deltaSnapshot`` = *[integer]* number of snapshots between reionization steps, i.e. computation of the ionization and photoionization rate fields; [vertical tree walking: ``endSnapshot`` - ``startSnapshot``, horizontal tree walking: 1]
- ``gridsize`` = *[integer]* gridsize of the density fields
- ``boxsize`` =  *[double]* length of the simulation box in h^-1 Mpc

**[Cosmology]**

- ``OM0`` = *[double]* matter density parameter
- ``OB0`` = *[double]* baryon density parameter
- ``OL0`` = *[double]* lambda density parameter
- ``HUBBLE_CONSTANT`` = *[double]* Hubble constant H = 100*h km/s/Mpc

**[StarFormation]**

- ``timestepModel`` = *[integer, optional]* flag to use either the default model with a constant fraction of gas turned into stars (as described in paper I) [value: 0], when comparing a simulation with a different time steps width but same spacing to a simulation using the default model [value: 1], or a model that assumes a constant star formation rate and is time step independent [value: 2]; if not defined, the default model is used
- ``timestepModel1_rescaleFactor`` = *[float, optional]* factor with which time steps of the default simulation needs to be multiplied with to obtain the time steps of the simulation that is compared to the default simulation
- ``timestepModel2_deltaTimeInMyr`` = *[float, optional]* time in Myr over which a fraction ``starFormationEfficiency`` of gas is turned into stars
- ``starFormationEfficiency`` = *[double]* maximum star formation efficiency [typical values: 0.01-0.03]

**[SNfeedback]**

- ``doDelayedSNfeedback`` = *[integer]* delayed [value: 1] or instantaneous [value: 0] supernova feedback
- ``SNenergyFractionIntoWinds`` = *[double]* supernova energy fraction that drives winds and causes gas ejection from galaxies [typical values: 0.05-0.3]

**[RadiativeFeedback]**

- ``doRadfeedback`` = *[integer]*  include [value: 1] or not include [value: 0] radiative feedback when computing the evolution of galaxy properties
- ``radfeedbackModel`` = radiative feedback model identifier; possible options are: MIN, SOBACCHI, TEMPEVOL, MJEANS
- ``ionThreshold`` = *[double]* ionization threshold above which a cell is considered as ionized [typical value: 0.5]
- ``tempIonGas`` = *[float]* temperature to which gas is heated upon ionization (Note for TEMPEVOL model: for M_c = M_F ``tempIonGas`` is a fourth of the temperature to which gas is heated upon ionization; for M_c = 8 M_F as indicated)
- ``muGas`` = *[float]* average particle mass in units of a proton mass [value: 0.59]

**[Reionization]**

- ``doReionization`` = *[integer]*
- ``cifogIniFile`` = path to ``cifogIniFile.ini``
- ``reionizationModel`` = flag to use either the self-consistent computed ionization field [flag: LOCAL] or impose the evolution found in Gnedin (2000) [GNEDIN]
- ``stellarPopulationSynthesisModel`` = stellar population synthesis model identifier which determines the number of ionizing photons; possible options are (suffix 'cont' indicates that star formation is assumed to be constinous across a timestep instead of being a delta function at the time of the snapshot): S99, S99cont, BPASS, BPASScont
- ``fescModel`` = escape fraction of ionizing photons model identifier; possible options are: CONST (constant fesc value defined under fescConst), MHDEC (fesc decreases with halo mass with boundary conditions defined under fescMH), MHINC  (fesc increases with halo mass with boundary conditions defined under fescMH), SN (fesc scales with the gas fraction ejected by supernovae feedback and is normalised by a factor which is given by ``fesc`` under fescConst)

**[fescConst]**

- ``fesc`` = *[double]* ionizing escape fraction value for CONST model, or normalisation factor for SN model

**[fescMH]**

- ``MHlow`` = *[double]* lowest halo mass where fesc is either 1 (MHDEC) or effectively 0 (MHINC)
- ``MHhigh`` = *[double]* highet halo mass where fesc is either 1 (MHINC) or effectively 0 (MHDEC)
- ``fescLow`` = *[double]* fesc value for the lowest halo mass
- ``fescHigh`` = *[double]* fesc value for the highest halo mass

**[Output]**

- ``type`` = *[integer]* flag to determine the format of the vertical outputs or tree files: number of trees, number of galaxies in tree 1, galaxies in tree 1, number of galaxies in tree 2, galaxies in tree 2, ... [value: 1]; number of trees, number of galaxies in tree 1, number of galaxies in tree 2, ..., galaxies in tree 1, galaxies in tree 2, ... [value: 2]
- ``horizontalOutput`` = *[integer]* write [value: 1] or do not write [value: 0] horizontal outputs, i.e. properties of all galaxies in a snapshot
- ``numSnapsToWrite`` = *[integer]* number of snapshots for which horizontal outputs should be written
- ``snapList`` = *[list of integers]* snapshot numbers for which horizontal outputs should be written [example: 12 25 34 38 42 46 51 54 56 58 62 64 69]
- ``verticalOutput`` = *[integer]* write [value: 1] or do not write [value: 0] vertical outputs or tree files constaining properties of galaxies in trees
- ``percentageOfTreesToWrite`` = *[integer]* percentage of trees to be written [default: 100]
- ``outputFile`` = path for directory where output files are to be written

Output Files
============

ASTRAEUS produces two types of output files:

- **horizontal:** file contains only galaxies at a chosen snapshot (redshift).
- **vertical:** file contains the local-horizontal merger trees with the computed galactic properties. 

Both output files are binary (and to date have not been adjusted for ) and have the following formats:

Horizontal
----------

Each galaxy has the following structure:

- *[float]* ``scalefactor``: scale factor 
- *[float]* ``pos[3]``: x, y, z coordinates in the simulation box in comoving Mpc/h
- *[float]* ``vel[3]``: velocities in x, z, z direction in physical km/s (peculiar)
- *[float]* ``Mvir``: halo mass in Msun/h
- *[float]* ``Mvir_prog``: sum of all progenitor halo masses in Msun/h
- *[float]* ``Rvir``: virial radius in comoving kpc/h
- *[float]* ``velDisp``: velocity dispersion in physical km/s
- *[float]* ``velMax``: maximum circular velocity in physical km/s
- *[float]* ``spin``: halo spin parameter
- *[float]* ``scalefactorLastMajorMerger``: scale factor of the last major merger (mass ratio > 0.3)
- *[float]* ``MgasIni``: initial gas mass (after gas accretion and radiative feedback, before star formation and supernoave feedback) in Msun/h
- *[float]* ``Mgas``: final gas mass (after star formation and supernovae feedback) in Msun/h
- *[float]* ``Mstar``: stellar mass in Msun/h
- *[float]* ``feff``: star formation efficiency 
- *[float]* ``fg``: fraction of the gas mass given by the cosmological ratio that a halo can retain in the presence of the UV background
- *[float]* ``zreion``: redshift when cell where galaxy is located became reionized
- *[float]* ``photHI_bg``: photoionization rate in s^-1 when cell where galaxy is located became reionized
- *[float]* ``stellarmasshistory[OUTPUTSNAPNUMBER]``: star formation history with each entry in the array listing the stellar mass form at the respective snapshot (redshift)

Vertical
--------

Each galaxy has the following structure:

- *[float]* ``scalefactor``: scale factor 
- *[float]* ``pos[3]``: x, y, z coordinates in the simulation box in comoving Mpc/h
- *[float]* ``vel[3]``: velocities in x, z, z direction in physical km/s (peculiar)
- *[float]* ``Mvir``: halo mass in Msun/h
- *[float]* ``Mvir_prog``: sum of all progenitor halo masses in Msun/h
- *[float]* ``Rvir``: virial radius in comoving kpc/h
- *[float]* ``velDisp``: velocity dispersion in physical km/s
- *[float]* ``velMax``: maximum circular velocity in physical km/s
- *[float]* ``spin``: halo spin parameter
- *[float]* ``scalefactorLastMajorMerger``: scale factor of the last major merger (mass ratio > 0.3)
- *[float]* ``MgasIni``: initial gas mass (after gas accretion and radiative feedback, before star formation and supernoave feedback) in Msun/h
- *[float]* ``Mgas``: final gas mass (after star formation and supernovae feedback) in Msun/h
- *[float]* ``Mstar``: stellar mass in Msun/h
- *[float]* ``feff``: star formation efficiency 
- *[float]* ``fg``: fraction of the gas mass given by the cosmological ratio that a halo can retain in the presence of the UV background
- *[float]* ``zreion``: redshift when cell where galaxy is located became reionized
- *[float]* ``photHI_bg``: photoionization rate in s^-1 when cell where galaxy is located became reionized
  
An example for a reading routines for the output files can be found in the analysis program linked below.
  
Analysis
========

The tree outputs generated with ``ASTRAEUS`` can be analysed using our analysis code `here <https://github.com/annehutter/astraeus/tree/master/analysis>`__.
