================
NEMO 4.0 Release
================

.. contents::

----------
What's new
----------

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

+-------------------+--------------------------------------------------------+------------------------------------+
| Name              | Purpose                                                | References                         |
+===================+=================+======================================+====================================+
| ``CANAL``         | East-west periodic canal of variable size with several |                                    |
|                   | initial states and associated geostrophic currents     |                                    |
|                   | (zonal jets or vortex).                                |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+
| ``ICEDYN``        | East-west + north-south periodic channel.              |                                    |
|                   | The common configuration includes an AGRIF zoom (1:3)  |                                    |
|                   | in the middle of the basin to test how an ice patch is |                                    |
|                   | advected through it but one can also test the          |                                    |
|                   | advection schemes (Prather and Ultimate-Macho) by      |                                    |
|                   | removing the ``key_agrif`` in the CPP keys.            |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+
| ``ISOMIP``        | Simple box configuration with an iceshelf with simple  | `Hunter 2006`_                     |
|                   | geometry on top.                                       |                                    |
|                   | The purpose of this test case is to evaluate the       |                                    |
|                   | impact of various schemes and new development with     |                                    |
|                   | iceshelf cavities.                                     |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+
| ``LOCK_EXCHANGE`` | Classical fluid dynamics experiment that has been      | - `Haidvogel and Beckmann 1999`_   |
|                   | adapted for testing advection schemes in ocean         | - `Burchard and Bolding 2002`_     |
|                   | circulation models.                                    | - `Ilıcak 2012`_                   |
|                   | This experiment can in particular illustrate the       |                                    |
|                   | impact of different choices of numerical schemes       |                                    |
|                   | and/or subgrid closures on spurious interior mixing.   |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+
| ``OVERFLOW``      | Adapted from the non-rotating overflow configuration   | - `Haidvogel and Beckmann 1999`_   |
|                   | Illustrates the impact of different choices of         | - `Ilıcak 2012`_                   |
|                   | numerical schemes and/or subgrid closures on spurious  |                                    |
|                   | interior mixing close to bottom topography.            |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+
| ``VORTEX``        | Illustrates the propagation of an anticyclonic eddy    | - `Debreu 2012`_                   |
|                   | over a Beta plan and flat bottom.                      | - `Penven 2006`_                   |
|                   | It is implemented here with an online refined          | - `Spall and Holland 1991`_        |
|                   | subdomain (thanks to AGRIF library) out of which the   |                                    |
|                   | vortex propagates and serves as a benchmark to         |                                    |
|                   | diagnose nesting errors.                               |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+
| ``WAD``           | Set of simple closed basin geometries for testing the  |                                    |
|                   | wetting and drying capabilities.                       |                                    |
|                   | Examples range from a closed channel with EW linear    |                                    |
|                   | bottom slope to a parabolic EW channel with a Gaussian |                                    |
|                   | ridge.                                                 |                                    |
+-------------------+--------------------------------------------------------+------------------------------------+

-----------
Improvments
-----------

Core components
===============

Passive tracer TOP and biogeochemical PISCES components
-------------------------------------------------------

- The passive tracers transport component was redesigned toward a modular structure and
  users can enable each module directly through logical flags in namelist_top (no more Fortran macros!).
- TOP on-line user documentation is available on NEMO Trac platform (TOP-UserQuickGuide_)
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
.. _TOP-UserQuickGuide:          http://forge.ipsl.jussieu.fr/nemo/wiki/WorkingGroups/top-dg/TOP-UserQuickGuide

.. _Hunter 2006:                 http://staff.acecrc.org.au/~bkgalton/ISOMIP/test_cavities.pdf
.. _Brodeau 2017:                http://doi.org/10.1175/JPO-D-16-0169.1
.. _Haidvogel and Beckmann 1999: http://hdl.handle.net/10013/epic.11761
.. _Burchard and Bolding 2002:   http://www.researchgate.net/publication/258128069_GETM_A_General_Estuarine_Transport_Model_Scientific_Documentation
.. _Ilıcak 2012:                 http://doi.org/10.1016/j.ocemod.2011.10.003
.. _Debreu 2012:                 http://doi.org/10.1016/j.ocemod.2012.03.003
.. _Penven 2006:                 http://doi.org/10.1016/j.ocemod.2005.05.002
.. _Spall and Holland 1991:      http://www.researchgate.net/publication/232101325_A_Nested_Primitive_Equation_Model_for_Oceanic_Applications
.. _Holland 2012:                http://doi.org/10.1175/JCLI-D-11-00078.1
.. _Lupkes 2012:                 http://doi.org/10.1029/2012JD017630
.. _Lupkes 2014:                 http://doi.org/10.1002/2014JD022418
.. _Kimmritz 2016:               http://doi.org/10.1016/j.ocemod.2016.03.004
.. _Orr and Epitalon 2015:       http://doi.org/10.5194/gmd-8-485-2015
