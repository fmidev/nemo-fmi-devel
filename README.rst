NEMO_ for **Nucleus for European Modelling of the Ocean** is a state-of-the-art modelling framework for research activities and forecasting services in ocean and climate sciences,
developed in a sustainable way by a European consortium since 2008.

Overview
========

The NEMO_ ocean model has 3 major components:

- **OCE** models the ocean [thermo]dynamics and solves the primitive equations (``./src/OCE``)
- **SI3** simulates seaice [thermo]dynamics, brine inclusions and subgrid-scale thickness variations (``./src/ICE``)
- **TOP** models the [on|off]line oceanic tracers transport and biogeochemical processes (``./src/TOP``)

These physical core engines are described in their respective `reference documentation`_ that must be cited for any work related to their use.

Assets and ready-to-use solutions
=================================

Not only does the NEMO_ framework model the ocean circulation, it offers various features to enable

- Create embedded zoom seamlessly thanks to 2-way nesting package AGRIF_.
- Opportunity to integrate an `external biogeochemistry model`_ 
- Versatile data_assimilation (``./src/OBS``)
- Generation of diagnostics through effective XIOS_ system
- Roll-out Earth system modeling with coupling interface based on OASIS_

Several built-in configurations are provided to evaluate the skills and performances of the model which can be used as templates for setting up a new configurations (``./cfgs``).

NEMO user can also checkout available idealized test cases that address specific physical processes(``./tests``).

A set of utilities is also provided to {pre,post}process your data (``./tools``).

Documentation and references
============================

A walkthrough tutorial illustrates how to get code dependencies, compile and execute NEMO (``./INSTALL.rst``) . 

Reference manuals and quick start guide can be build from source and exported to HTML or PDF formats (``./doc``) or downloaded directly from the NEMO_ website.

=========== =================== ===================
 Component    Reference Manual   Quick start
=========== =================== ===================
 OCE          `NEMO manual`_     `NEMO guide`_
 SI3          `SI3 manual`_
 TOP          `TOP manual`_      
=========== =================== ===================

Since 2014 the project has a `Special Issue`_ in the open-access journal Geoscientific Model Development (GMD) from the European Geosciences Union (EGU).
The main scope is to collect relevant manuscripts covering various topics and to provide a single portal to assess the model potential and evolution.

Used by a wide audience, numerous `associated projects`_ have been carried out and extensive `bibliography`_ published.

NEMO Consortium
===============

The NEMO Consortium pulling together 5 European institutes (CMCC_, CNRS_, MOI_, `Met Office`_ and NERC_) plans the sustainable development in order to keep a reliable evolving framework since 2008.

It defines the NEMO strategy that is implemented by the System Team on a yearly basis in order to release a new version almost every four years.

When the need arises, `working groups`_ are created or resumed to gather the community expertise for advising on the development activities.

..  _external biogeochemistry model : http://forge.ipsl.jussieu.fr/nemo/wiki/WorkingGroups/TOP/TOP-UserQuickGuide
.. _working groups : https://forge.ipsl.jussieu.fr/nemo/wiki/WorkingGroups
.. _reference documentation : https://www.nemo-ocean.eu/bibliography/documentation/
.. _Special Issue : http://www.geosci-model-dev.net/special_issue40.html
.. _associated projects : https://www.nemo-ocean.eu/projects/
.. _bibliography : https://www.nemo-ocean.eu/wp-content/plugins/wp-bibtexbrowser/bibtexbrowser.php?bib=nemo.bib
