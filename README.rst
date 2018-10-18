.. include:: .global.rst

:Release: |release|
:Date:    |today|
            
NEMO for *Nucleus for European Modelling of the Ocean* is a state-of-the-art modelling framework of 
ocean related engines for oceanographic research, operational oceanography, seasonal forecast and 
[paleo]climate studies.

.. contents::
	:local:

Overview
========

Distributed under CeCILL license (GNU GPL compatible - see ``./LICENSE``),
the NEMO ocean model has 3 major components:

- :rstblue:`OPA` is fundamental to all users by modelling the ocean [thermo]dynamics and 
  solving the primitive equations         (``./src/OCE``); :cite:`madec_bk08`
- :rstgrey:`SI`\ :rstgreysup:`3` for sea-ice simulates ice [thermo]dynamics, brine inclusions and 
  subgrid-scale thickness variations      (``./src/ICE``); :cite:`gmd-8-2991-2015,vancoppenolle200933`
- :rstgreen:`TOP-PISCES` models biogeochemistry with TOP for the on/offline oceanic tracers transport and 
  PISCES for the biogeochemical processes (``./src/MBG``). :cite:`gmd-8-2465-2015`

These physical engines are described in their respective reference publications that must be cited for
any work related to their use.

Applications and capabilities
=============================

Not only does the NEMO framework model the ocean circulation,
it offers various features to enable

- 2-way nesting package `AGRIF`_ to create :doc:`embedded zooms <zooms>` seamlessly
- Flexible biogeochemistry with :doc:`online coarsening <coarsening>` and
  opportunity to integrate an :doc:`alternative model <tracers>`
- Versatile :doc:`data assimilation interface <data_assimilation>`

lation a efficient XIOS_ server for outputing diagnostics a coupled via OASIS_ to alternative components or other models to enable Earth system modelling.

| Several :doc:`builtins configurations <reference_configurations>` are provided to assess the skills and
	performances of the model which can be used as templates for setting up a new configuration (``./cfgs``).
| The end user could also find some :doc:`idealised test cases <test_cases>` on the web to serve as examples and
	to study particular processes (``./tests``).

A set of tools is also provided to setup your own configuration and [pre|post]process your data (``./tools``).

Literature
==========

| The NEMO reference manual and a quick start guide can be generated from Latex and RST source code files
	(``./doc``), either in PDF or in HTML format, but it might require some additional installations.
| In any case, both formats are available online: `HTML`_ | `PDF`_

| Since 2014 the project has a `Special Issue`_ in the open-access journal Geoscientific Model Development (GMD) 
	from the European Geosciences Union (EGU).
	The main scope is to collect relevant manuscripts which cover a wide variety of topics like
	process studies, new parameterizations, implementation of new model features and new NEMO configurations.
| Also it provides a single portal to search, discover and understand about
	the NEMO modelling framework potential and evolution and to submit their contributions. 

Community development
=====================

| The NEMO Consortium gathering 6 European institutes (`CMCC`_, `CNRS`_, `MOI`_, `Met Office`_ and `NERC`_)
	organises the sustainable development in order to keep a reliable evolving system since 2008.
| It defines the multi-year development strategy which is implemented by the NEMO System Team.

When the need arises, `Working Groups`_ are created or resumed to gather the expertise in the community in order to
focus the development work on a specific subject or major component of the framework.

How to cite NEMO
================

..	bibliography:: references.bib
	:all:
	:style: unsrt

.. _HTML:           http://www.nemo-ocean.eu/doc
.. _PDF:            http://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf
