=====================
Build a configuration
=====================

.. include:: .global.rst

.. contents::
	:local:
	:depth: 1
      
.. role:: underline 
   :class: underline

Official configurations
=======================

| NEMO is distributed with some reference configurations allowing both the user to set up a first application and
  the developer to validate their developments.
| :underline:`The NEMO System Team is in charge of these configurations`.

+----------------------+-----+-----+-----+--------+-------+-------------------------------+
|                      | OPA | SI3 | TOP | PISCES | AGRIF | Inputs                        |
+======================+=====+=====+=====+========+=======+===============================+
| `AGRIF_DEMO`_        |  X  |  X  |     |        |   X   | - `AGRIF_DEMO_v4.0.tar`_      |
|                      |     |     |     |        |       | - `ORCA2_ICE_v4.0.tar`_       |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `AMM12`_             |  X  |     |     |        |       | `AMM12_v4.0.tar`_             |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `C1D_PAPA`_          |  X  |     |     |        |       | `INPUTS_C1D_PAPA_v4.0.tar`_   |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `GYRE_BFM`_          |  X  |     |  X  |        |       | ``-``                         |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `GYRE_PISCES`_       |  X  |     |  X  |   X    |       | ``-``                         |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `ORCA2_ICE_PISCES`_  |  X  |  X  |  X  |   X    |       | - `ORCA2_ICE_v4.0.tar`_       |
|                      |     |     |     |        |       | - `INPUTS_PISCES_v4.0.tar`_   |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `ORCA2_OFF_PISCES`_  |     |     |  X  |   X    |       | - `INPUTS_PISCES_v4.0.tar`_   |
|                      |     |     |     |        |       | - `ORCA2_OFF_v4.0.tar`_       |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `ORCA2_OFF_TRC`_     |     |     |  X  |        |       | `ORCA2_OFF_v4.0.tar`_         |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `ORCA2_SAS_ICE`_     |     |  X  |     |        |       | - `ORCA2_ICE_v4.0.tar`_       |
|                      |     |     |     |        |       | - `INPUTS_SAS_v4.0.tar`_      |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+
| `SPITZ12`_           |  X  |  X  |     |        |       | `SPITZ12_v4.0.tar`_           |
+----------------------+-----+-----+-----+--------+-------+-------------------------------+

----------
AGRIF_DEMO
----------

.. image:: _static/AGRIF_DEMO.jpg

``AGRIF_DEMO`` is based on the ``ORCA2_LIM3_PISCES`` global 2° configuration but
it includes 3 online nested grids that demonstrate the overall capabilities of AGRIF in a realistic context,
including nesting sea ice models.

The configuration includes a 1:1 grid in the Pacific and two successively nested grids with odd and
even refinement ratios over the Arctic ocean.
The finest grid spanning the whole Svalbard archipelago is of particular interest to check that
sea ice coupling is done properly.
The 1:1 grid, used alone, is used as a benchmark to check that the solution is not corrupted by grid exchanges.

Note that since grids interact only at the baroclinic time level,
numerically exact results can not be achieved in the 1:1 case.
One has to switch to a fully explicit in place of a split explicit free surface scheme in order to
retrieve perfect reproducibility.

Corresponding ``AGRIF_FixedGrids.in`` file is given by::

	2
	42 82 49 91 1 1 1
	122 153 110 143 4 4 4
	0
	1
	38 80 71 111 3 3 3
	0

-----
AMM12
-----

``AMM12`` for *Atlantic Margin Model 12kms* is a `regional model`_ covering the Northwest European Shelf domain on
a regular lat-lon grid at approximately 12km horizontal resolution.
The key ``key_amm_12km`` is used to create the correct dimensions of the AMM domain.

| This configuration tests several features of NEMO functionality specific to the shelf seas.
| In particular, the AMM uses s-coordinates in the vertical rather than z-coordinates and is forced with
  tidal lateral boundary conditions using a flather boundary condition from the BDY module (``key_bdy``).

The AMM configuration uses the GLS (``key_zdfgls``) turbulence scheme,
the VVL non-linear free surface (``key_vvl``) and time-splitting (``key_dynspg_ts``).

In addition to the tidal boundary condition, the model may also take open boundary conditions from
a North Atlantic model.
Boundaries may be completely ommited by removing the BDY key (key_bdy) in ``./cfgs/AMM12/cpp_AMM12_fcm``.

Sample surface fluxes, river forcing and a sample initial restart file are included to test a realistic model run.
The Baltic boundary is included within the river input file and is specified as a river source.
Unlike ordinary river points the Baltic inputs also include salinity and temperature data.

--------
C1D_PAPA
--------

``C1D_PAPA`` is a 1D configuration (one water column called NEMO1D, activated with CPP key ``key_c1d``),
located at the `PAPA station 145W-50N <http://www.pmel.noaa.gov/OCS/Papa/index-Papa.shtml>`_.

| NEMO1D is useful to test vertical physics in NEMO
  (turbulent closure scheme, solar penetration, interaction ocean/atmosphere.,...)
| Size of the horizontal domain is 3x3 grid points.

This reference configuration uses a 75 vertical levels grid (1m at the surface),
the GLS (key_zdfgls) turbulence scheme with K-epsilon closure and the CORE BULK formulae.
The atmospheric forcing comes from ECMWF operational analysis with a modification of the long and short waves flux.
This set has been rescaled at a frequency of 1h. 1 year is simulated in outputs,
see below (June,15 2010 to June,14 2011)

`Reffray 2015`_ describes some tests on vertical physic using this configuration.

The inputs tar file includes:

- forcing files covering the years 2010 and 2011 (``forcing_PAPASTATION_1h_y201*.nc``)
- initialization file for June,15 2010 deduced from observed data and Levitus 2009 climatology
  (``init_PAPASTATION_m06d15.nc``)
- surface chlorophyll file (``chlorophyll_PAPASTATION.nc``) deduced from Seawifs data.

--------
GYRE_BFM
--------

``GYRE_BFM`` is the same configuration as `GYRE_PISCES`_, except that PISCES is replaced by
BFM biogeochemichal model in coupled mode.

-----------
GYRE_PISCES
-----------

| Idealized configuration representing double gyres in the North hemisphere, Beta-plane with
  a regular grid spacing at 1° horizontal resolution (and possible use as a benchmark by
  easily inscreasing grid size), 101 vertical levels, forced with analytical heat, freshwater and
  wind-stress fields.
| This configuration is coupled to `PISCES biogeochemical model`_.

Running GYRE as a benchmark
---------------------------

This simple configuration can be used as a benchmark since it is easy to increase resolution
(and in this case no physical meaning of outputs):

1. Choose the grid size

   In ``./cfgs/GYRE/EXP00``, edit your ``namelist_cfg`` file to change the ``jp_cfg``, ``jpi``, ``jpj``,
   ``jpk`` variables in &namcfg:

	+------------+---------+---------+---------+------------------+---------------+
	| ``jp_cfg`` | ``jpi`` | ``jpj`` | ``jpk`` | Number of points | Equivalent to |
	+============+=========+=========+=========+==================+===============+
	| 1          | 30      | 20      | 101     | 60600            | GYRE 1°       |
	+------------+---------+---------+---------+------------------+---------------+
	| 25         | 750     | 500     | 101     | 37875000         | ORCA 1/2°     |
	+------------+---------+---------+---------+------------------+---------------+
	| 50         | 1500    | 1000    | 101     | 151500000        | ORCA 1/4°     |
	+------------+---------+---------+---------+------------------+---------------+
	| 150        | 4500    | 3000    | 101     | 1363500000       | ORCA 1/12°    |
	+------------+---------+---------+---------+------------------+---------------+
	| 200        | 6000    | 4000    | 101     | 2424000000       | ORCA 1/16°    |
	+------------+---------+---------+---------+------------------+---------------+

2. In `namelist_cfg` again, avoid problems in the physics (and results will not be meaningful in terms of physics) by setting `nn_bench = 1` in &namctl

.. code-block:: fortran
   
   nn_bench    =    1     !  Bench mode (1/0): CAUTION use zero except for bench

3. If you increase domain size, you may need to decrease time-step (for stability) by changing `rn_rdt` value in &namdom (i.e. for `jp_cfg = 150`, ORCA12 equivalent, use `rn_rdt = 1200`)

.. code-block:: fortran
   
   rn_rdt      = 1200.     !  time step for the dynamics

4. Optional, in order to increase the number of MPI communication for benchmark purposes:
   you can change the number of sub-timesteps computed in the time-splitting scheme each iteration.
   First change the list of active CPP keys for your experiment,
   in `cfgs/"your configuration name"/cpp_"your configuration name".fcm`:
   replace ``key_dynspg_flt by key_dynspg_ts`` and recompile/create your executable again
   
   .. code-block:: fortran
   
   makenemo [...] add_key 'key_dynspg_ts' del_key 'key_dynspg_flt'

In your ``namelist_cfg`` file, edit the &namsplit namelist by adding the following line: 

.. code-block:: fortran
   
   nn_baro       =    30               !  Number of iterations of barotropic mode/

``nn_baro = 30`` is a kind of minimum (we usually use 30 to 60).
So than increasing the ``nn_baro`` value will increase the number of MPI communications.

The GYRE CPP keys, namelists and scripts can be explored in the ``GYRE`` configuration directory
(``./cfgs/GYRE`` and ``./cfgs/GYRE/EXP00``).

Find monthly mean outputs of 1 year run here:
http://prodn.idris.fr/thredds/catalog/ipsl_public/reee451/NEMO_OUT/GYRE/catalog.html

----------------
ORCA2_ICE_PISCES
----------------

ORCA is the generic name given to global ocean configurations.
Its specificity lies on the horizontal curvilinear mesh used to overcome the North Pole singularity found for
geographical meshes.
SI3 (Sea Ice Integrated Initiative) is a thermodynamic-dynamic sea ice model specifically designed for
climate studies.
A brief description of the model is given here.

:underline:`Space-time domain`

The horizontal resolution available through the standard configuration is ORCA2.
It is based on a 2 degrees Mercator mesh, (i.e. variation of meridian scale factor as cosinus of the latitude).
In the northern hemisphere the mesh has two poles so that the ratio of anisotropy is nearly one everywhere.
The mean grid spacing is about 2/3 of the nominal value: for example it is 1.3 degrees for ORCA2.
Other resolutions (ORCA4, ORCA05 and ORCA025) are running or under development within specific projects.
In the coarse resolution version (i.e. ORCA2 and ORCA4) the meridional grid spacing is increased near
the equator to improve the equatorial dynamics.
Figures in pdf format of mesh and bathymetry can be found and downloaded here.
The sea-ice model runs on the same grid.

The vertical domain spreads from the surface to a depth of 5000m.
There are 31 levels, with 10 levels in the top 100m.
The vertical mesh is deduced from a mathematical function of z ([[AttachmentNum(1)]]).
The ocean surface corresponds to the w-level k=1, and the ocean bottom to the w-level k=31.
The last T-level (k=31) is thus always in the ground.The depths of the vertical levels and
the associated scale factors can be viewed.
Higher vertical resolution is used in ORCA025 and ORCA12 (see `DRAKKAR project <http://www.drakkar-ocean.eu>`_).

The time step depends on the resolution. It is 1h36' for ORCA2 so that there is 15 time steps in one day.

:underline:`Ocean Physics (for ORCA2)`

- horizontal diffusion on momentum: the eddy viscosity coefficient depends on the geographical position.
  It is taken as 40000 $m^2/s$, reduced in the equator regions (2000 $m^2/s$) excepted near the western boundaries.
- isopycnal diffusion on tracers: the diffusion acts along the isopycnal surfaces (neutral surface) with
  a eddy diffusivity coefficient of 2000 $m^2/s$.
- Eddy induced velocity parametrization with a coefficient that depends on the growth rate of
  baroclinic instabilities (it usually varies from 15 $m^2/s$ to 3000 $m^2/s$).
- lateral boundary conditions : zero fluxes of heat and salt and no-slip conditions are applied through
  lateral solid boundaries.
- bottom boundary condition : zero fluxes of heat and salt are applied through the ocean bottom.
  The Beckmann [19XX] simple bottom boundary layer parameterization is applied along continental slopes.
  A linear friction is applied on momentum.
- convection: the vertical eddy viscosity and diffusivity coefficients are increased to 1 $m^2/s$ in case of
  static instability.
- forcings: the ocean receives heat, freshwater, and momentum fluxes from the atmosphere and/or the sea-ice.
  The solar radiation penetrates the top meters of the ocean.
  The downward irradiance I(z) is formulated with two extinction coefficients [Paulson and Simpson, 1977],
  whose values correspond to a Type I water in Jerlov's classification (i.e the most transparent water)

ORCA2_ICE_PISCES is a reference configuration with the following characteristics:

- global ocean configuration
- based on a tri-polar ORCA grid, with a 2° horizontal resolution
- 31 vertical levels
- forced with climatological surface fields
- coupled to the sea-ice model SI3.
- coupled to TOP passive tracer transport module and `PISCES biogeochemical model`_.

:underline:`AGRIF demonstrator`

| From the ``ORCA2_ICE_PISCES`` configuration, a demonstrator using AGRIF nesting can be activated.
  It includes the global ``ORCA2_ICE_PISCES`` configuration and a nested grid in the Agulhas region.
| To set up this configuration, after extracting NEMO:

- Build your AGRIF configuration directory from ORCA2_ICE_PISCES, with the key_agrif CPP key activated:

.. code-block:: console
                
	$ ./makenemo -r 'ORCA2_ICE_PISCES' -n 'AGRIF' add_key 'key_agrif'

- Using the ``ORCA2_ICE_PISCES`` input files and namelist, AGRIF test configuration is ready to run

:underline:`On-The-Fly Interpolation`

| NEMO allows to use the interpolation on the fly option allowing to interpolate input data during the run.
  If you want to use this option you need files giving informations on weights, which have been created.
| You can find at http://prodn.idris.fr/thredds/catalog/ipsl_public/reee512/ORCA2_ONTHEFLY/WEIGHTS/catalog.html
  2 weights files `bil_weights` for scalar field (bilinear interpolation) and `bic_weights` for
  vector field (bicubic interpolation).
| The data files used are `COREII forcing <http://data1.gdfl.noaa.gov/nomads/forms/mom4/COREv2>`_ extrapolated on
  continents, ready to be used for on the fly option:
  `COREII`_ forcing files extrapolated on continents

----------------
ORCA2_OFF_PISCES
----------------

``ORCA2_OFF_PISCES`` uses the ORCA2 configuration in which the `PISCES biogeochemical model`_ has been activated in
standalone using the dynamical fields that are pre calculated.

See `ORCA2_ICE_PISCES`_ for general description of ORCA2.

The input files for PISCES are needed, in addition the dynamical fields are used as input.
They are coming from a 2000 years of an ORCA2_LIM climatological run using ERA40 atmospheric forcing.

-------------
ORCA2_OFF_TRC
-------------

``ORCA2_OFF_TRC`` uses the ORCA2_LIM configuration in which the tracer passive transport module TOP has been
activated in standalone using the dynamical fields that are pre calculated.

See `ORCA2_ICE_PISCES`_ for general description of ORCA2.

In ``namelist_top_cfg``, different passive tracers can be activated ( cfc11, cfc12, sf6, c14, age ) or my-trc,
a user-defined tracer.

The dynamical fields are used as input, they are coming from a 2000 years of an ORCA2_LIM climatological run using
ERA40 atmospheric forcing.

-------------
ORCA2_SAS_ICE
-------------

``ORCA2_SAS_ICE`` is a demonstrator of the SAS ( Stand-alone Surface module ) based on ORCA2_LIM configuration.

The standalone surface module allows surface elements such as sea-ice, iceberg drift and surface fluxes to
be run using prescribed model state fields.
For example, it can be used to inter-compare different bulk formulae or adjust the parameters of
a given bulk formula

See `ORCA2_ICE_PISCES`_ for general description of ORCA2.

Same input files as `ORCA2_ICE_PISCES`_ are needed plus fields from a previous ORCA2_LIM run.

More informations on input and configuration files in `NEMO Reference manual`_.

-------
SPITZ12
-------

``SPITZ12``

Unsupported configurations
==========================

Other configurations are developed and used by some projects with "NEMO inside",
these projects are welcome to publicize it here: http://www.nemo-ocean.eu/projects/add-project

:underline:`Obviously these "projects configurations" are not under the NEMO System Team's responsibility`.

.. _regional model:               http://www.tandfonline.com/doi/pdf/10.1080/1755876X.2012.11020128
.. _AMM12_v4.0.tar:               http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/AMM12_v4.0.tar
.. _PISCES biogeochemical model:  http://www.geosci-model-dev.net/8/2465/2015
.. _INPUTS_PISCES_v4.0.tar:       http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/INPUTS_PISCES_v4.0.tar
.. _ORCA2_OFF_v4.0.tar:           http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/ORCA2_OFF_v4.0.tar
.. _ORCA2_ICE_v4.0.tar:           http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/ORCA2_ICE_v4.0.tar
.. _INPUTS_SAS_v4.0.tar:          http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/INPUTS_SAS_v4.0.tar
.. _NEMO Reference manual:        http://forge.ipsl.jussieu.fr/nemo/doxygen/index.html?doc=NEMO
.. _INPUTS_C1D_PAPA_v4.0.tar:     http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/INPUTS_C1D_PAPA_v4.0.tar
.. _Reffray 2015:                 http://www.geosci-model-dev.net/8/69/2015
.. _COREII:                       http://prodn.idris.fr/thredds/catalog/ipsl_public/reee512/ORCA2_ONTHEFLY/FILLED_FILES/catalog.html
.. _SPITZ12_v4.0.tar:             http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/SPITZ12_v4.0.tar
.. _AGRIF_DEMO_v4.0.tar:          http://prodn.idris.fr/thredds/fileServer/ipsl_public/romr005/Online_forcing_archives/AGRIF_DEMO_v4.0.tar
