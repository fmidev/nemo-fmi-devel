======================
What's new in NEMO 4.0
======================

.. include:: .global.rst

.. contents::
	:local:
      
Original sea-ice component SI\ :sup:`3`\
========================================

| Meeting in 2017, the `sea ice working group`_ decided to gather sea ice developers and
  scientific expertise within the NEMO community in order to build a unified model.
  The replacement of LIM (Louvain-la-Neuve Ice Model) has required the backport of
  the desired functionalities from each of CICE, GELATO and LIM into a common code in the NEMO framework.
| This new model has been named SI\ :sup:`3`\ ("sea ice cubed") and
  means "**Sea Ice modelling Integrated Initiative**".

- Physics

  * Landfast ice: simple grounding with a "release" stress (not fully operational)
  * Lateral melting
  * Melt ponds: constant or `Holland 2012`_ formulation (and soon topographic melt ponds)
  * Ice-atm. drags from `Lupkes 2012`_ (depending on ice concentration) or `Lupkes 2014`_
    (Depending on sea ice concentration and atm. stability)
- Numerics

  * Advection: Ultimate-Macho scheme
  * Rheology: adaptive EVP (`Kimmritz 2016`_)
  * Coupling interface: conductivity as surface forcing instead of heat fluxes (Met Office requirement)
- Performance

  * All thermodynamics in 1D
  * Reduced mpp communications
- Users & developers friendly

  * Comprehensive set of outputs (universal units and understandable names + includes limp)
  * New architecture and namelist
  * All processes can be decoupled from each other (switch on/off)
  * Ice categories bounds can be defined by the user or set automatically
  * For open boundaries, the number of ice categories from the forcing model can be different
    from the number of categories in the regional simulation
  * Fully compatible with AGRIF

First Test Cases
================

Define and install a separate repository for test cases to all easy contributions from the NEMO Users Community

Core components
===============

Passive tracer TOP and biogeochemical PISCES components
-------------------------------------------------------

- The passive tracers transport component was redesigned toward a modular structure and
  users can enable each module directly through logical flags in namelist_top (no more Fortran macros!).
- TOP on-line user documentation is available on NEMO Trac platform (`TOP User Quick Guide`_)
- TOP currently accounts for the following 5 modules:

  * ``CFC`` contains inorganic carbon tracers (CFC11/CFC12/SF6)
  * ``MY_TRC`` is a template for new modules (or external couplings)
  * ``AGE`` deals with water age tracking
  * ``C14`` as a radiocarbon passive tracer
  * ``PISCES`` the companion ecosystem model
- A generalized infrastructure was developed to handle the prescription of either surface, coastal, or
  open boundaries conditions for each passive tracer.
- PISCES model contains new developments and modifications:

  * Particulate Organic Carbon (POC) component comes with a new liability scheme,
    while the former Kriest parameterisation was superseded;
  * A complex iron chemistry scheme is now available, with an improved description of ligands for
    the marine iron cycle
  * Carbonate chemistry is based on MOCSY 2.0 routines (see `Orr and Epitalon 2015`_),
    by complying also with CMIP6 standards.
  * Ecosystem components can be optionally modelled by means of explicit nutrient quotas (PISCES-QUOTA)

AGRIF (embedded zooms)
----------------------

The NEMO 4.0 includes new capabilities, configurations and test cases with AGRIF:

.. role:: underline
	:class: underline

:underline:`New capabilities from NEMO 3.6 to NEMO 4.0`

AGRIF is continuously maintained so that it could be activated with all NEMO components (OPA, sea-ice, TOP).
Depending on NEMO version, it is nevertheless not the case so that some options may not be compatible with
the use of online grid refinement.
Check out the table below to know the status according to the NEMO release you may use.

:underline:`Status of available options with AGRIF (if not listed, option is compatible with AGRIF)`:

+--------------------------------------------------------+----------------+---------------------+
|                                                        | NEMO 3.6       | NEMO 4.0            |
+========================================================+================+=====================+
| LIM2                                                   | yes            | ``-``               |
+--------------------------------------------------------+----------------+---------------------+
| LIM3/SI3                                               | no             | yes                 |
+--------------------------------------------------------+----------------+---------------------+
| TOP                                                    | yes            | yes                 |
+--------------------------------------------------------+----------------+---------------------+
| GLS vertical mixing                                    | no             | yes                 |
+--------------------------------------------------------+----------------+---------------------+
| z*                                                     | no             | yes                 |
+--------------------------------------------------------+----------------+---------------------+
| z~                                                     | no             | no                  |
+--------------------------------------------------------+----------------+---------------------+
| Lagrangian icebergs                                    | no             | no                  |
+--------------------------------------------------------+----------------+---------------------+
| East-west periodic and/or north fold bcs in zooms      | no             | no                  |
+--------------------------------------------------------+----------------+---------------------+
| Online timing                                          | no             | no                  |
+--------------------------------------------------------+----------------+---------------------+
| Stochastic parameterization                            | no             | no                  |
+--------------------------------------------------------+----------------+---------------------+
| Vertical coordinate change in zooms (``key_vertical``) | no             | yes, but not tested |
+--------------------------------------------------------+----------------+---------------------+
| Number of ghost cells                                  | 1 (hard coded) | 3 (parameter)       |
+--------------------------------------------------------+----------------+---------------------+

[Important notice concerning the change of ghost cells number]

The default number of ghost cells (i.e. the number of cells that serve as open boundary data provision) has been
increased from 1 to 3 in NEMO 4.0.
This allows to properly handle boundary conditions for numerical schemes that
have a discretization order greater than 2.
On the user point of view this does not change anything++ except in the definition of level 1 grids in
the ``AGRIF_FixedGrids.in`` file.
In order to retrieve exactly the position of a nested grid in NEMO 4.0 one has to shift indices by
2 points to the south-west.
Taking the ``ICEDYN`` example above for NEMO 4.0, the "old" NEMO 3.6 corresponding file would contain::

	1
	36 65 36 65 3 3 3
	0

++ Child grid output files are now greater by 4 points in each direction.

- Now compatible with new sea ice component and z* coordinate
- Extended ghost cells area to properly handle scheme with spatial order >2
- Added vertical refinement (beta)
- Nesting tools for setup now up to date and working

Adding of 3 Reference Configurations
------------------------------------

- ``AGRIF_DEMO``: 2 interlocked zooms (1:4 & 1:3) in the Nordic Seas + 1 zoom (1:1) at the equator
- ``ORCA2_OFF_TRC``: a benchmark simulation environment to deal with inert carbon tracers dynamics by
  exploiting the offline coupling with NEMO.
- ``SPITZ12``: regional configuration around the Svalbard archipelago.

Misc. improvments
=================

Physics
-------

- Bulk formulae : move to aerobulk package (`Brodeau 2017`_), i.e. NCAR, COARE and ECMWF bulk
  (remove Clio and MFS bulk)
- Fix for tracer conservation with split explicit free surface
- Wetting and drying
- iso-neutral mixing (iso and triad operators): add the Method of Stabilizing Correction (MSC)
  (more accurate calculation) + add a bilaplacian case
- Lateral physics (LDF): scale aware setting of eddy viscosity and diffusivity
- Wave coupling: large scale wave interaction process added in momentum and tracer equations
- Remove the acceleration of convergence

Numerics
--------

- Added tidal self attraction and loading either read from a file or from usual "scalar" approximation
- Vertical physics (ZDF) (modularity, share shear production calculation between TKE and GKS,
  removal of all ZDF CPP keys, removal of avmu & avmv, minimization of MPP comm.: ~15 removed)
- Remove the split-explicit ZDF scheme for both TRA and DYN
- Lateral physics (LDF): simplification of user interface and removal of CPP keys
- Add a 4th order centered (CEN) and Flux Corrected Transport (FCT) tracer advection
  (using a 4th compact in the vertical)
- Generalised lbc_lnk and lbc_nfd
- Configuration interface completely rewritten (DOM module mainly suppressed,
  and in place: domain_cfg.nc file, or usr_def module)
- Vorticity: 2 new energy conserving scheme:  ENT with Coriolis defined at T-point
  (better for Flux form) and EET a variant of EEN where e3t is used instead of e3f
  (solved the issue with e3f specification but is not enstrophy conserving)
- Wave coupling: coupled interface to external wave model

Performances
------------

- MPI Message passing recoded to reduce number of MPI communications (suppression of redundant communications,
  gather multiple communications into one)
- Back to standard dynamical allocation (remove of wrk_alloc/dealloc statements)
- XIOS software for IOs version 2 as default, and optionally available for restarts
- Unify mppini
- Use non uniform jpi/jpj with dynamic allocation to avoid ghost rows/columns

Environment
-----------

- Revised structure of namelist_ref/_cfg and default reference values.
- Reorganisation of SVN repository to be compliant with usual directory tree and facilitate building of
  NEMO executable
- Improvements of reliability through automatic and regular testing of the changes made in repository

.. _sea ice working group:       http://forge.ipsl.jussieu.fr/nemo/wiki/WorkingGroups/SI3
.. _TOP User Quick Guide:        http://forge.ipsl.jussieu.fr/nemo/wiki/WorkingGroups/top-dg/TOP-UserQuickGuide

.. _Brodeau 2017:                http://doi.org/10.1175/JPO-D-16-0169.1
.. _Holland 2012:                http://doi.org/10.1175/JCLI-D-11-00078.1
.. _Lupkes 2012:                 http://doi.org/10.1029/2012JD017630
.. _Lupkes 2014:                 http://doi.org/10.1002/2014JD022418
.. _Kimmritz 2016:               http://doi.org/10.1016/j.ocemod.2016.03.004
.. _Orr and Epitalon 2015:       http://doi.org/10.5194/gmd-8-485-2015
