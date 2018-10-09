===================
Interfacing options
===================

.. contents::
	:local:
	:depth: 1
           
Embedded zooms with AGRIF
=========================

.. contents::
   :local:

--------
Overview
--------

AGRIF (Adaptive Grid Refinement In Fortran) is a library that allows the seamless space and time refinement over
rectangular regions in NEMO.
Refinement factors can be odd or even (usually lower than 5 to maintain stability).
Interaction between grid is "two-ways" in the sense that the parent grid feeds the child grid open boundaries and
the child grid provides volume averages of prognostic variables once a given number of time step is completed.
These pages provide guidelines how to use AGRIF in NEMO.
For a more technical description of the library itself, please refer to http://agrif.imag.fr.

-----------
Compilation
-----------

Activating AGRIF requires to append the cpp key ``key_agrif`` at compilation time: 

.. code-block:: sh

	./makenemo add_key 'key_agrif'

Although this is transparent to users, the way the code is processed during compilation is different from
the standard case:
a preprocessing stage (the so called "conv" program) translates the actual code so that
saved arrays may be switched in memory space from one domain to an other.

--------------------------------
Definition of the grid hierarchy
--------------------------------

An additional text file ``AGRIF_FixedGrids.in`` is required at run time.
This is where the grid hierarchy is defined.
An example of such a file, here taken from the ``ICEDYN`` test case, is given below::

	1
	34 63 34 63 3 3 3
	0

The first line indicates the number of zooms (1).
The second line contains the starting and ending indices in both directions on the root grid
(imin=34 imax=63 jmin=34 jmax=63) followed by the space and time refinement factors (3 3 3).
The last line is the number of child grid nested in the refined region (0).
A more complex example with telescoping grids can be found below and
in the ``AGRIF_DEMO`` reference configuration directory.

[Add some plots here with grid staggering and positioning ?]

When creating the nested domain, one must keep in mind that the child domain is shifted toward north-east and
depends on the number of ghost cells as illustrated by the (attempted) drawing below for nbghostcells=1 and
nbghostcells=3.
The grid refinement is 3 and nxfin is the number of child grid points in i-direction.  

.. image:: _static/agrif_grid_position.jpg

Note that rectangular regions must be defined so that they are connected to a single parent grid.
Hence, defining overlapping grids with the same refinement ratio will not work properly,
boundary data exchange and update being only performed between root and child grids.
Use of east-west periodic or north-fold boundary conditions is not allowed in child grids either.
Defining for instance a circumpolar zoom in a global model is therefore not possible. 

-------------
Preprocessing
-------------

Knowing the refinement factors and area, a ``NESTING`` pre-processing tool may help to create needed input files
(mesh file, restart, climatological and forcing files).
The key is to ensure volume matching near the child grid interface,
a step done by invoking the ``Agrif_create_bathy.exe`` program.
You may use the namelists provided in the ``NESTING`` directory as a guide.
These correspond to the namelists used to create ``AGRIF_DEMO`` inputs.

----------------
Namelist options
----------------

Each child grid expects to read its own namelist so that different numerical choices can be made
(these should be stored in the form ``1_namelist_cfg``, ``2_namelist_cfg``, etc... according to their rank in
the grid hierarchy).
Consistent time steps and number of steps with the chosen time refinement have to be provided.
Specific to AGRIF is the following block:

.. code-block:: fortran

	!-----------------------------------------------------------------------
	&namagrif      !  AGRIF zoom                                            ("key_agrif")
	!-----------------------------------------------------------------------
	   ln_spc_dyn    = .true.  !  use 0 as special value for dynamics
	   rn_sponge_tra = 2880.   !  coefficient for tracer   sponge layer [m2/s]
	   rn_sponge_dyn = 2880.   !  coefficient for dynamics sponge layer [m2/s]
	   ln_chk_bathy  = .false. !  =T  check the parent bathymetry
	/             

where sponge layer coefficients have to be chosen according to the child grid mesh size.
The sponge area is hard coded in NEMO and applies on the following grid points:
2 x refinement factor (from i=1+nbghostcells+1 to i=1+nbghostcells+sponge_area) 

----------   
References
----------

`Debreu, L., P. Marchesiello, P. Penven and G. Cambon, 2012: Two-way nesting in split-explicit ocean models: Algorithms, implementation and validation. Ocean Modelling, 49-50, 1-21. <http://doi.org/10.1016/j.ocemod.2012.03.003>`_

`Penven, P., L. Debreu, P. Marchesiello and J. C. Mc Williams, 2006: Evaluation and application of the ROMS 1-way embedding procedure to the central california upwelling system. Ocean Modelling, 12, 157-187. <http://doi.org/10.1016/j.ocemod.2005.05.002>`_

`Spall, M. A. and W. R. Holland, 1991: A Nested Primitive Equation Model for Oceanic Applications. J. Phys. Ocean., 21, 205-220. <https://doi.org/10.1175/1520-0485(1991)021\<0205:ANPEMF\>2.0.CO;2>`_

----

On line biogeochemistry coarsening
==================================

.. contents::
   :local:

.. role:: underline 
   :class: underline

------------
Presentation
------------

A capacity of coarsening physics to force a BGC model coupled to NEMO has been developed.
This capacity allow to run 'online' a BGC model coupled to OCE-SI3 with a lower resolution,
to reduce the CPU cost of the BGC model, while preserving the effective resolution of the dynamics.

A presentation is available [attachment:crs_wiki_1.1.pdf​ here], where the methodology is presented.

-----------------------------------------------------
What is available and working for now in this version
-----------------------------------------------------

[To be completed]

----------------------------------------------
Description of the successful validation tests
----------------------------------------------

[To be completed]

------------------------------------------------------------------
What is not working yet with on line coarsening of biogeochemistry
------------------------------------------------------------------

[To be completed]

''should include precise explanation on MPI decomposition problems too''

---------------------------------------------
How to set up and use on line biogeochemistry
---------------------------------------------

:underline:`How to activate coarsening?`

To activate the coarsening, ``key_crs`` should be added to list of CPP keys.
This key will only activate the coarsening of dynamics.

Some parameters are available in the namelist_cfg:

.. code-block:: fortran

	               !   passive tracer coarsened online simulations
	!-----------------------------------------------------------------------
	   nn_factx    = 3         !  Reduction factor of x-direction
	   nn_facty    = 3         !  Reduction factor of y-direction
	   nn_msh_crs  = 0         !  create (=1) a mesh file or not (=0)
	   nn_crs_kz   = 3         ! 0, volume-weighted MEAN of KZ
	                           ! 1, MAX of KZ
	                           ! 2, MIN of KZ
	                           ! 3, 10^(MEAN(LOG(KZ)) 
	                           ! 4, MEDIANE of KZ 
	   ln_crs_wn   = .false.   ! wn coarsened (T) or computed using horizontal divergence ( F )
	                           !                           !
	   ln_crs_top = .true.     !coarsening online for the bio
	/

- Only ``nn_factx = 3`` is available and the coarsening only works for grids with a T-pivot point for
  the north-fold lateral boundary condition (ORCA025, ORCA12, ORCA36, ...).
- ``nn_msh_crs = 1`` will activate the generation of the coarsened grid meshmask.
- ``nn_crs_kz`` is the operator to coarsen the vertical mixing coefficient. 
- ``ln_crs_wn``

  - when ``key_vvl`` is activated, this logical has no effect;
    the coarsened vertical velocities are computed using horizontal divergence.
  - when ``key_vvl`` is not activated,

    - coarsened vertical velocities are computed using horizontal divergence (``ln_crs_wn = .false.``) 
    - or coarsened vertical velocities are computed with an average operator (``ln_crs_wn = .true.``)
- ``ln_crs_top = .true.``: should be activated to run BCG model in coarsened space;
  so only works when ``key_top`` is in the cpp list and eventually ``key_pisces`` or ``key_my_trc``.

:underline:`Choice of operator to coarsene KZ`

A sensiblity test has been done with an Age tracer to compare the different operators.
The 3 and 4 options seems to provide the best results.

Some results can be found [xxx here]

:underline:`Example of xml files to output coarsened variables with XIOS`

In the [attachment:iodef.xml iodef.xml]  file, a "nemo" context is defined and
some variable defined in [attachment:file_def.xml file_def.xml] are writted on the ocean-dynamic grid.  
To write variables on the coarsened grid, and in particular the passive tracers,
a "nemo_crs" context should be defined in [attachment:iodef.xml iodef.xml] and
the associated variable are listed in [attachment:file_crs_def.xml file_crs_def.xml ].

:underline:`Passive tracers tracers initial conditions`

When initial conditions are provided in NetCDF files, the field might be:

- on the coarsened grid
- or they can be on another grid and
  interpolated `on-the-fly <http://forge.ipsl.jussieu.fr/nemo/wiki/Users/SetupNewConfiguration/Weight-creator>`_.
  Example of namelist for PISCES :
  
	.. code-block:: fortran

		!-----------------------------------------------------------------------
		&namtrc_dta      !    Initialisation from data input file
		!-----------------------------------------------------------------------
		!
		   sn_trcdta(1)  = 'DIC_REG1'        ,        -12        ,  'DIC'     ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(2)  = 'ALK_REG1'        ,        -12        ,  'ALK'     ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(3)  = 'O2_REG1'         ,        -1         ,  'O2'      ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(5)  = 'PO4_REG1'        ,        -1         ,  'PO4'     ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(7)  = 'Si_REG1'         ,        -1         ,  'Si'      ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
	   	sn_trcdta(10) = 'DOC_REG1'        ,        -12        ,  'DOC'     ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(14) = 'Fe_REG1'         ,        -12        ,  'Fe'      ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(23) = 'NO3_REG1'        ,        -1         ,  'NO3'     ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
	   	rn_trfac(1)   =   1.0e-06  !  multiplicative factor
		   rn_trfac(2)   =   1.0e-06  !  -      -      -     -
		   rn_trfac(3)   =  44.6e-06  !  -      -      -     -
	   	rn_trfac(5)   = 122.0e-06  !  -      -      -     -
		   rn_trfac(7)   =   1.0e-06  !  -      -      -     -
		   rn_trfac(10)  =   1.0e-06  !  -      -      -     -
	   	rn_trfac(14)  =   1.0e-06  !  -      -      -     -
		   rn_trfac(23)  =   7.6e-06  !  -      -      -     -
		
	   	cn_dir        =  './'      !  root directory for the location of the data files

:underline:`PISCES forcing files`

They might be on the coarsened grid.

:underline:`Perspectives`

For the future, a few options are on the table to implement coarsening for biogeochemistry in 4.0 and
future releases.
Those will be discussed in Autumn 2018

----

Coupling with other models (OASIS, SAS, ...)
============================================

NEMO currently exploits OASIS-3-MCT to implement a generalised coupled interface
(`Coupled Formulation <http://forge.ipsl.jussieu.fr/nemo/doxygen/node50.html?doc=NEMO>`_).
It can be used to interface with most of the European atmospheric GCM (ARPEGE, ECHAM, ECMWF, Ha- dAM, HadGAM, LMDz),
as well as to WRF (Weather Research and Forecasting Model), and to implement the coupling of
two independent NEMO components, ocean on one hand and sea-ice plus other surface processes on the other hand
(`Standalone Surface Module - SAS <http://forge.ipsl.jussieu.fr/nemo/doxygen/node46.html?doc=NEMO>`_).

To enable the OASIS interface the required compilation key is ``key_oasis3``.
The parameters to set are in sections ``namsbc_cpl`` and in case of using of SAS also in section ``namsbc_sas``.

----

With data assimilation
======================

.. contents::
   :local:

The assimilation interface to NEMO is split into three modules.
- OBS for the observation operator
- ASM for the application of increments and model bias correction (based on the assimilation increments).
- TAM the tangent linear and adjoint model.

Please see the `NEMO reference manual`_ for more details including information about the input file formats and
the namelist settings.

--------------------------------------
Observation and model comparison (OBS)
--------------------------------------

The observation and model comparison code (OBS) reads in observation files (profile temperature and salinity,
sea surface temperature, sea level anomaly, sea ice concentration, and velocity) and
calculates an interpolated model equivalent value at the observation location and nearest model timestep.
The resulting data are saved in a feedback file (or files).
The code was originally developed for use with the NEMOVAR data assimilation code, but
can be used for validation or verification of model or any other data assimilation system.
This is all controlled by the namelist.
To build with the OBS code active ``key_diaobs`` must be set. 

More details in the `NEMO reference manual`_ chapter 12.

Standalone observation operator (SAO)
-------------------------------------

The OBS code can also be run after a model run using saved NEMO model data.
This is accomplished using the standalone observation operator (SAO)
(previously known the offline observation operator).

To build the SAO use makenemo.
This means compiling NEMO once (in the normal way) for the chosen configuration.
Then include ``SAO`` at the end of the relevant line in ``cfg.txt`` file.
Then recompile with the replacement main program in ``./src/SAO``.
This is a special version of ``nemogcm.F90`` (which doesn't run the model, but reads in the model fields, and
observations and runs the OBS code.
See section 12.4 of the `NEMO reference manual`_.

-----------------------------------
Apply assimilation increments (ASM)
-----------------------------------

The ASM code adds the functionality to apply increments to the model variables:
temperature, salinity, sea surface height, velocity and sea ice concentration.
These are read into the model from a NetCDF file which may be produced by separate data assimilation code.
The code can also output model background fields which are used as an input to data assimilation code.
This is all controlled by the namelist nam_asminc.
To build the ASM code ``key asminc`` must be set.

More details in the `NEMO reference manual`_ chapter 13.

--------------------------------
Tangent linear and adjoint (TAM)
--------------------------------

This is the tangent linear and adjoint code of NEMO which is useful to 4D VAR assimilation.

----

Inputs-Outputs (using XIOS)
===========================

.. contents::
   :local:

| Output of diagnostics in NEMO is usually done using XIOS.
  This is an efficient way of writing diagnostics because the time averaging, file writing and even
  some simple arithmetic or regridding is carried out in parallel to the NEMO model run.
| This page gives a basic introduction to using XIOS with NEMO.
  Much more information is available from the XIOS homepage above and from the `NEMO reference manual`_.

Use of XIOS for diagnostics is activated using the pre-compiler key ``key_iomput``.
The default version of XIOS is the 2.0 release. 

------------------------------
Extracting and installing XIOS
------------------------------

1. Install the NetCDF4 library.
   If you want to use single file output you will need to compile the HDF & NetCDF libraries to allow parallel IO.
2. Download the version of XIOS that you wish to use.
   The recommended version is now XIOS 2.0:
   
	.. code-block:: console

		$ svn co ​http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.0 xios-2.0

   and follow the instructions in `XIOS documentation`_ to compile it.
   If you find problems at this stage, support can be found by subscribing to the `XIOS users mailing list`_ and
   sending a mail message to it. 

---------
Namelists
---------

XIOS is controlled using xml input files that should be copied to your model run directory before
running the model.
The exact setup differs slightly between 1.0 and 2.0 releases.

An ``iodef.xml`` file is still required in the run directory.
For XIOS 2.0 the ``field_def.xml`` file has been further split into ``field_def-oce.xml`` (for physics),
``field_def-ice.xml`` (for ice) and ``field_def-bgc.xml`` (for biogeochemistry).
Also the definition of the output files has been moved from the ``iodef.xml`` file into
separate ``file_definition.xml`` files which are included in the ``iodef.xml`` file.
Note that the ``domain_def.xml`` file is also different for XIOS 2.0.

-----
Modes
-----

Detached Mode
-------------

In detached mode the XIOS executable is executed on separate cores from the NEMO model.
This is the recommended method for using XIOS for realistic model runs.
To use this mode set ``using_server`` to ``true`` at the bottom of the ``iodef.xml`` file:

.. code-block:: xml

	<variable id="using_server" type="boolean">true</variable>

Make sure there is a copy (or link to) your XIOS executable in the working directory and
in your job submission script allocate processors to XIOS.

Attached Mode
-------------

In attached mode XIOS runs on each of the cores used by NEMO.
This method is less efficient than the detached mode but can be more convenient for testing or
with small configurations.
To activate this mode simply set ``using_server`` to false in the ``iodef.xml`` file

.. code-block:: xml

	<variable id="using_server" type="boolean">false</variable>

and don't allocate any cores to XIOS.
Note that due to the different domain decompositions between XIOS and NEMO if
the total number of cores is larger than the number of grid points in the j direction then the model run will fail.

------------------------------
Adding new diagnostics to NEMO
------------------------------

If you want to add a NEMO diagnostic to the NEMO code you will need to do the following:

1. Add any necessary code to calculate you new diagnostic in NEMO
2. Send the field to XIOS using ``CALL iom_put( 'field_id', variable )`` where ``field_id`` is a unique id for
   your new diagnostics and variable is the fortran variable containing the data.
   This should be called at every model timestep regardless of how often you want to output the field.
   No time averaging should be done in the model code. 
3. If it is computationally expensive to calculate your new diagnostic you should also use "iom_use" to
   determine if it is requested in the current model run. For example,
   
	.. code-block:: fortran

		IF iom_use('field_id') THEN
		   !Some expensive computation
	   	!...
		   !...
	   	iom_put('field_id', variable)
		ENDIF

4. Add a variable definition to the ``field_def.xml`` (or ``field_def-???.xml``) file
5. Add the variable to the ``iodef.xml`` or ``file_definition.xml`` file.

.. _NEMO reference manual:   http://forge.ipsl.jussieu.fr/nemo/doxygen/index.html?doc=NEMO
.. _XIOS documentation:      http://forge.ipsl.jussieu.fr/ioserver/wiki/documentation
.. _XIOS users mailing list: http://forge.ipsl.jussieu.fr/mailman/listinfo.cgi/xios-users
