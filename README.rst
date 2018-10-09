.. role:: rstblue
.. role:: rstgreen
.. role:: rstgrey
.. role:: rstgreysup(sup)

NEMO for *Nucleus for European Modelling of the Ocean* is a state-of-the-art modelling framework of 
ocean related engines for oceanographic research, operational oceanography, seasonal forecast and 
[paleo]climate studies.

Overview
========

Distributed under CeCILL license (GNU GPL compatible - see ``./LICENSE``),
the NEMO ocean model has 3 major components:

- :rstblue:`OPA` is fundamental to all users by modelling the ocean [thermo]dynamics and 
  solving the primitive equations         (``./src/OCE``); 
- :rstgrey:`SI`\ :rstgreysup:`3` for sea-ice simulates ice [thermo]dynamics, brine inclusions and 
  subgrid-scale thickness variations      (``./src/ICE``); 
- :rstgreen:`TOP-PISCES` models biogeochemistry with TOP for the on/offline oceanic tracers transport and 
  PISCES for the biogeochemical processes (``./src/MBG``).

These physical engines are described in their respective reference publications that must be cited for 
any work related to their use.

Applications and capabilities
=============================

Not only does the NEMO framework model the ocean circulation,
it offers various features to enable

- 2-way nesting package `AGRIF`_ to create embedded zooms seamlessly
- Flexible biogeomchemistry with online coarsening and possible integration of a customized model 
- Versatile data assimilation interface with 3 different modules
  (tangent linear, observational operators and increments)

lation a efficient XIOS_ server for outputing diagnostics a coupled via OASIS_ to alternative components or other models to enable Earth system modelling.

| Several builtins configurations are provided to assess the skills and performances of the model which
	can be used as templates for setting up a new configuration (``./cfgs``).
| The end user could also find some idealised test cases on the web to serve as examples and
	to study particular processes (``./tests``).

A set of tools is also provided to setup your own configuration and [pre|post]process your data (``./tools``).

Literature
==========

| The NEMO reference manual and a quick start guide can be generated from Latex and RST source code files
	(``./doc``), either in PDF or in HTML format, but it might require some additional installations.
| In any case, both formats are available online: `HTML`_ | `PDF`_

| Since 2014 the project has a `Special Issue`_ in the Geoscientific Model Development (GMD) open-access journal
	from the European Geosciences Union (EGU).
	The main scope is to collect relevant manuscripts which cover a wide variety of topics like
   process studies, new parameterizations, implementation of new model features and new NEMO configurations.
| Also it provides a single portal to search, discover and understand about
	the NEMO modelling framework potential and evolution and to submit their contributions. 

Community development
=====================

| The NEMO Consortium gathering 6 European institutes organises the sustainable development in order to 
	keep a reliable evolving system since 2008.
| It defines the multi-year development strategy which is implemented by the NEMO System Team.

`Working groups`_ are regularly created or resumed to gather the expertise in the NEMO community in order to 
focus the development work on a specific subject or major component of NEMO.

Definitions
===========

AGRIF
	*Adaptive Grid Refinement In Fortran*, 
	package for the integration of full adaptive mesh refinement features within 
	an existing multidimensional finite difference model

SI\ :sup:`3`\ 
	*Sea Ice Integrated Initiative*, 
	unified sea ice model merging functionalities from CICE, GELATO and LIM into the NEMO framework

OASIS
	*Ocean Atmosphere Sea Ice Soil*, 
	coupling software to synchronise numerical codes representing different components of the climate system

PISCES
	*Pelagic Interactions Scheme for Carbon and Ecosystem Studies*, 
	biogeochemical model simulating marine ecosystems, cycles of carbon and the main nutrients

TAM
	*Tangent linear and Adjoint Model*, 
	tools to analyse and control the NEMO dynamical core for a wide range of applications such as
	sensitivity analysis, parameter estimation, vectors computation or data assimilation.

TOP
	*Tracers in Ocean Paradigm*, 
	on/off-line oceanic tracers transport and biogeochemistry models

XIOS
	*XML Input Output Server*, 
	library dedicated to input/output management of climate code

.. _AGRIF:          http://agrif.imag.fr
.. _HTML:           http://www.nemo-ocean.eu/doc
.. _NEMO:           http://www.nemo-ocean.eu
.. _OASIS:          http://verc.enes.org/oasis
.. _PDF:            http://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf
.. _Special Issue:  http://www.geosci-model-dev.net/special_issue40.html
.. _Working groups: http://forge.ipsl.jussieu.fr/nemo/wiki/WorkingGroups
.. _XIOS:           http://forge.ipsl.jussieu.fr/ioserver
