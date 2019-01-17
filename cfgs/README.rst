REFERENCE CONFIGURATIONS
========================

NEMO is distributed with a set of reference configurations allowing both the user to set up his own first applications and the developer to test/validate his NEMO developments (using SETTE package).
*The NEMO System Team is in charge only for these configurations.*

Configurations developed by external research projects or initiatives that make use of NEMO are welcome to be publicized also through the NEMO website by filling up the form `Add project <http://www.nemo-ocean.eu/projects/add>`_

Available Configurations 
------------------------
====================== ===== ===== ===== ======== ======= ================================================
 Configuration                     Component(s)                            Input & Forcing File(s)
---------------------- ---------------------------------- ------------------------------------------------
 Name                   OPA   SI3   TOP   PISCES   AGRIF
====================== ===== ===== ===== ======== ======= ================================================
 AGRIF_DEMO_             X     X                     X     AGRIF_DEMO_v4.0.tar_, ORCA2_ICE_v4.0.tar_
 AMM12_                  X                                 AMM12_v4.0.tar_
 C1D_PAPA_               X                                 INPUTS_C1D_PAPA_v4.0.tar_
 GYRE_BFM_               X           X                     *none*
 GYRE_PISCES_            X           X      X              *none*
 ORCA2_ICE_PISCES_       X     X     X      X              ORCA2_ICE_v4.0.tar_, INPUTS_PISCES_v4.0.tar_
 ORCA2_OFF_PISCES_                   X      X              ORCA2_OFF_v4.0.tar_, INPUTS_PISCES_v4.0.tar_
 ORCA2_OFF_TRC_                      X                     ORCA2_OFF_v4.0.tar_
 ORCA2_SAS_ICE_                X                           ORCA2_ICE_v4.0.tar_, INPUTS_SAS_v4.0.tar_
 SPITZ12_                X     X                           SPITZ12_v4.0.tar_
====================== ===== ===== ===== ======== ======= ================================================

**How to compile an experiment from a reference configuration**

A user who wants to compile the ORCA2_ICE_PISCES_ reference configuration using makenemo should use the following, by selecting among available architecture file or providing a user defined one:


.. code-block:: console
                
        $ ./makenemo -r 'ORCA2_ICE_PISCES' -m 'my-fortran.fcm' -j '4'

A new EXP00 folder will be created within the selected reference configurations, namely ``trunk/cfgs/ORCA2_ICE_PISCES/EXP00``, where it will be necessary to uncompress the Input & Forcing Files listed in the above table.

Then it will be possible to launch the execution of the model through a runscript (opportunely adapted to the user system).

AGRIF_DEMO
----------

AGRIF_DEMO is based on the ORCA2_ICE_PISCES_ global configuration at 2° of resolution with the inclusion of 3 online nested grids to demonstrate the overall capabilities of AGRIF (Adaptive Grid Refinement In Fortran) in a realistic context (including the nesting of sea ice models).

The configuration includes a 1:1 grid in the Pacific and two successively nested grids with odd and even refinement ratios over the Arctic ocean, with the finest grid spanning the whole Svalbard archipelago that is of particular interest to test sea ice coupling.

The 1:1 grid can be used alone as a benchmark to check that the model solution is not corrupted by grid exchanges. 
Note that since grids interact only at the baroclinic time level, numerically exact results can not be achieved in the 1:1 case. Perfect reproducibility is obtained only by switching to a fully explicit setip instead of a split explicit free surface scheme.

AMM12
-----

AMM12 stands for *Atlantic Margin Model at 12 km* that is a regional configuration covering the Northwest European Shelf domain on a regular horizontal grid of ~12 km of resolution (see `O'Dea et al., 2012 <http://www.tandfonline.com/doi/pdf/10.1080/1755876X.2012.11020128>`_).

This configuration allows to tests several features of NEMO specifically addressed to the shelf seas. 
In particular, AMM12  accounts for vertical s-coordinates system, GLS turbulence scheme, tidal lateral boundary conditions using a flather scheme (see more in BDY).
Boundaries may be completely omitted by setting ``ln_bdy = .false.`` in ``nambdy``.

Sample surface fluxes, river forcing and an initial restart file are included to test a realistic model run (AMM12_v4.0.tar_).

Note that, the Baltic boundary is included within the river input file and is specified as a river source, but unlike ordinary river points the Baltic inputs also include salinity and temperature data.

C1D_PAPA
--------

C1D_PAPA is a 1D configuration for the `PAPA station <http://www.pmel.noaa.gov/OCS/Papa/index-Papa.shtml>`_ located in the northern-eastern Pacific Ocean at 50.1°N, 144.9°W. See `Reffray et al. (2015) <http://www.geosci-model-dev.net/8/69/2015>`_ for the description of its physical and numerical turbulent-mixing behaviour.

The water column setup, called NEMO1D, is activated with the inclusion of the CPP key ``key_c1d`` and has a horizontal domain of 3x3 grid points.

This reference configuration uses 75 vertical levels grid (1m at the surface), GLS turbulence scheme with K-epsilon closure and the NCAR bulk formulae.
Data provided with INPUTS_C1D_PAPA_v4.0.tar_ file account for :

- ``forcing_PAPASTATION_1h_y201[0-1].nc`` : ECMWF operational analysis atmospheric forcing rescaled to 1h (with long and short waves flux correction) for years 2010 and 2011
- ``init_PAPASTATION_m06d15.nc`` : Initial Conditions from observed data and Levitus 2009 climatology
- ``chlorophyll_PAPASTATION.nc`` : surface chlorophyll file from Seawifs data


GYRE_BFM
--------

GYRE_BFM shares the same physical setup of GYRE_PISCES_, but NEMO is coupled with the `BFM <http://www.bfm-community.eu/>`_ biogeochemical model as described in :trac:`source:/NEMO/trunk/cfgs/GYRE_BFM/README`.


GYRE_PISCES
-----------

GYRE_PISCES is an idealized configuration representing a Northern hemisphere double gyres system,  in the Beta-plane approximation with a regular 1° horizontal resolution and 31 vertical levels, which is coupled with `PISCES biogeochemical model`_. Analytical forcing for heat, freshwater and wind-stress fields are applied.  

This configuration act also as demonstrator of the **USER DEFINED setup** (``ln_read_cfg = .false.``) and grid setting are handled through the ``&namusr_def`` controls in namelist_cfg:

.. code-block:: fortran

  !-----------------------------------------------------------------------
  &namusr_def    !   GYRE user defined namelist
  !-----------------------------------------------------------------------
     nn_GYRE     =     1     !  GYRE resolution [1/degrees]
     ln_bench    = .false.   !  ! =T benchmark with gyre: the gridsize is kept constant
     jpkglo      =    31     !  number of model levels
  /

Note that, the default grid size is 30x20 grid points (with ``nn_GYRE = 1``) and vertical levels are set by ``jpkglo``. The specific code changes can be inspected at :trac:`source:/NEMO/trunk/src/OCE/USR` 

**Running GYRE as a benchmark** :  this simple configuration can be used as a benchmark since it is easy to increase resolution, with the drawback of getting results that have a very limited physical meaning.

GYRE grid resolution can be increased at runtime by setting a different value of ``nn_GYRE`` (integer multiplier scaling factor),  as described in the following table: 

=========== ========= ========== ============ ===================
``nn_GYRE``  *jpiglo*  *jpjglo*   ``jpkglo``   **Equivalent to**
=========== ========= ========== ============ ===================
 1           30        20         31           GYRE 1°
 25          750       500        101          ORCA 1/2°
 50          1500      1000       101          ORCA 1/4°
 150         4500      3000       101          ORCA 1/12°
 200         6000      4000       101          ORCA 1/16°
=========== ========= ========== ============ ===================

Note that,  it is necessary to set ``ln_bench = .true.`` in ``namusr_def`` to avoid problems in the physics computation and that the model timestep should be adequately rescaled. 

For example if ``nn_GYRE = 150``, equivalent to an ORCA 1/12° grid, the timestep should be set to 1200 seconds

.. code-block:: fortran
   
   rn_rdt      = 1200.     !  time step for the dynamics

Differently from previous versions of NEMO, the code uses by default  the time-splitting scheme and internally computes the number of sub-steps. 


ORCA2_ICE_PISCES
----------------

ORCA2_ICE_PISCES is a reference configuration for the global ocean with a 2°x2° curvilinear horizontal mesh and 31 vertical levels, distributed using z-coordinate system and with 10 levels in the top 100m.
ORCA is the generic name given to global ocean Mercator mesh, (i.e. variation of meridian scale factor as cosinus of the latitude), with two poles in the northern hemisphere so that the ratio of anisotropy is nearly one everywhere

In this configuration, the ocean dynamical core  is coupled to  

- **ICE**, namely SI3 (Sea Ice Integrated Initiative) a thermodynamic-dynamic sea ice model specifically designed for climate studies.
- **TOP**, passive tracer transport module and `PISCES biogeochemical model`_

All components share the same grid.

The model is forced with CORE-II normal year atmospheric forcing and it uses the NCAR bulk formulae.

**Ocean Physics configuration**

- *horizontal diffusion on momentum*: the eddy viscosity coefficient depends on the geographical position. It is taken as 40000 m^2/s, reduced in the equator regions (2000 m^2/s) excepted near the western boundaries.
- *isopycnal diffusion on tracers*: the diffusion acts along the isopycnal surfaces (neutral surface) with an eddy diffusivity coefficient of 2000 m^2/s.
- *Eddy induced velocity parametrization* with a coefficient that depends on the growth rate of baroclinic instabilities (it usually varies from 15 m^2/s to 3000 m^2/s).
- *lateral boundary conditions* : zero fluxes of heat and salt and no-slip conditions are applied through lateral solid boundaries.
- *bottom boundary condition* : zero fluxes of heat and salt are applied through the ocean bottom.
  The Beckmann [19XX] simple bottom boundary layer parameterization is applied along continental slopes.
  A linear friction is applied on momentum.
- *convection*: the vertical eddy viscosity and diffusivity coefficients are increased to 1 m^2/s in case of static instability.
- *time step* is 5760sec (1h36') so that there is 15 time steps in one day.



**AGRIF demonstrator**

From the ORCA2_ICE_PISCES configuration, a demonstrator using AGRIF nesting can be activated that includes a nested grid in the Agulhas region.

To set up this configuration, after extracting NEMO:

Build your AGRIF configuration directory from ORCA2_ICE_PISCES, with the key_agrif CPP key activated:

.. code-block:: console
                
        $ ./makenemo -r 'ORCA2_ICE_PISCES' -n 'AGRIF' add_key 'key_agrif'

By using the input files and namelists for ORCA2_ICE_PISCES, the AGRIF test configuration is ready to run.


ORCA2_OFF_PISCES
----------------

ORCA2_OFF_PISCES  shares the same general offline configuration of ORCA2_ICE_TRC, but only PISCES model is an active component of TOP.


ORCA2_OFF_TRC
-------------

ORCA2_OFF_TRC is based on the ORCA2 global ocean configuration (see `ORCA2_ICE_PISCES`_ for general description) along with the tracer passive transport module (TOP), but dynamical fields are pre-calculated and read with specific time frequency.

This enables for an offline coupling of TOP components, here specifically inorganic carbon compounds (cfc11, cfc12, sf6, c14) and water age module (age). See ``namelist_top_cfg`` to inspect the selection of each component with the dedicated logical keys.

Pre-calculated dynamical fields are provided to NEMO using the namelist ``&namdta_dyn``  in ``namelist_cfg``, in this case with a 5 days frequency (120 hours):

.. code-block:: fortran

  !-----------------------------------------------------------------------
  &namdta_dyn    !   offline ocean input files                            (OFF_SRC only)
  !-----------------------------------------------------------------------
     ln_dynrnf       =  .false.    !  runoffs option enabled (T) or not (F)
     ln_dynrnf_depth =  .false.    !  runoffs is spread in vertical (T) or not (F)
     cn_dir      = './'      !  root directory for the ocean data location
     !___________!_________________________!___________________!___________!_____________!________!___________!__________________!__________!_______________!
     !           !  file name              ! frequency (hours) ! variable  ! time interp.!  clim  ! 'yearly'/ ! weights filename ! rotation ! land/sea mask !
     !           !                         !  (if <0  months)  !   name    !   (logical) !  (T/F) ! 'monthly' !                  ! pairing  !    filename   !
     sn_tem      = 'dyna_grid_T'           ,       120         , 'votemper'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_sal      = 'dyna_grid_T'           ,       120         , 'vosaline'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_mld      = 'dyna_grid_T'           ,       120         , 'somixhgt'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_emp      = 'dyna_grid_T'           ,       120         , 'sowaflup'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_fmf      = 'dyna_grid_T'           ,       120         , 'iowaflup'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_ice      = 'dyna_grid_T'           ,       120         , 'soicecov'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_qsr      = 'dyna_grid_T'           ,       120         , 'soshfldo'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_wnd      = 'dyna_grid_T'           ,       120         , 'sowindsp'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_uwd      = 'dyna_grid_U'           ,       120         , 'uocetr_eff',  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_vwd      = 'dyna_grid_V'           ,       120         , 'vocetr_eff',  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_wwd      = 'dyna_grid_W'           ,       120         , 'wocetr_eff',  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_avt      = 'dyna_grid_W'           ,       120         , 'voddmavs'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_ubl      = 'dyna_grid_U'           ,       120         , 'sobblcox'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_vbl      = 'dyna_grid_V'           ,       120         , 'sobblcoy'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
  /

Input dynamical fields for this configuration (ORCA2_OFF_v4.0.tar_) comes from a 2000 years long climatological simulation of ORCA2_ICE using ERA40 atmospheric forcing.

Note that, this configuration default uses linear free surface (``ln_linssh = .true.``) assuming that model mesh is not varying in time and it includes the bottom boundary layer parameterization (``ln_trabbl = .true.``) that requires the provision of bbl coefficients through ``sn_ubl`` and ``sn_vbl`` fields.

It is also possible to activate PISCES model (see ORCA2_OFF_PISCES_) or a user defined set of tracers and source-sink terms with ``ln_my_trc = .true.`` (and adaptation of :trac:`source:/NEMO/trunk/src/TOP/MY_TRC` routines).

In addition, the offline module (OFF) allows for the provision of further fields:

1. **River runoff** can be provided to TOP components by setting ``ln_dynrnf = .true.`` and by including an input datastream similarly to the following:

.. code-block:: fortran

     sn_rnf      = 'dyna_grid_T'           ,       120         , 'sorunoff'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''

2. **VVL dynamical fields**, in the case input data were produced by a dyamical core using variable volume (``ln_linssh = .false.``) it necessary to provide also diverce and E-P at before timestep by including input datastreams similarly to the following

.. code-block:: fortran

     sn_div       = 'dyna_grid_T'           ,       120         ,    'e3t'     ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''
     sn_empb      = 'dyna_grid_T'           ,       120         , 'sowaflupb'  ,  .true.   , .true. , 'yearly'  , ''               , ''       , ''


More details can be found by inspecting the offline data manager at :trac:`source:/NEMO/trunk/src/OFF/dtadyn.F90`


ORCA2_SAS_ICE
-------------

ORCA2_SAS_ICE is a demonstrator of the Stand-Alone Surface (SAS) module and it relies on ORCA2 global ocean configuration (see `ORCA2_ICE_PISCES`_ for general description).

The standalone surface module allows surface elements such as sea-ice, iceberg drift, and surface fluxes to be run using prescribed model state fields.
It can profitably be used to compare different bulk formulae or adjust the parameters of a given bulk formula.

More informations about SAS can be found in NEMO manual.

SPITZ12
-------

SPITZ12 is a regional configuration around the Svalbard archipelago at 1/12° of horizontal resolution and 75 vertical levels. See `Rousset et al. (2015) <https://www.geosci-model-dev.net/8/2991/2015/>`_ for more details.

This configuration references to year 2002, with atmospheric forcing provided every 2 hours using NCAR bulk formulae, while lateral boundary conditions for dynamical fields have 3 days time frequency.


.. _AGRIF_DEMO_v4.0.tar:          http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/AGRIF_DEMO_v4.0.tar
.. _AMM12_v4.0.tar:               http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/AMM12_v4.0.tar
.. _PISCES biogeochemical model:  http://www.geosci-model-dev.net/8/2465/2015
.. _INPUTS_PISCES_v4.0.tar:       http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/INPUTS_PISCES_v4.0.tar
.. _ORCA2_OFF_v4.0.tar:           http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/ORCA2_OFF_v4.0.tar
.. _ORCA2_ICE_v4.0.tar:           http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/ORCA2_ICE_v4.0.tar
.. _INPUTS_SAS_v4.0.tar:          http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/INPUTS_SAS_v4.0.tar
.. _INPUTS_C1D_PAPA_v4.0.tar:     http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/INPUTS_C1D_PAPA_v4.0.tar
.. _SPITZ12_v4.0.tar:             http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/SPITZ12_v4.0.tar

.. _COREII:                       http://prodn.idris.fr/thredds/catalog/ipsl_public/reee512/ORCA2_ONTHEFLY/FILLED_FILES/catalog.html

