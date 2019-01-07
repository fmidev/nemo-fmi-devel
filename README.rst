:Release: |release|
:Date:    |today|

`NEMO`_ for *Nucleus for European Modelling of the Ocean* is a state-of-the-art modelling framework,
for research activities and forecasting services in ocean and climate sciences,
developed in a sustainable way by a European consortium since 2008.

.. contents::
	:local:

Overview
========

The NEMO ocean model has 3 major components:

- |OPA| models the ocean [thermo]dynamics and solves the primitive equations
  (``./src/OCE``) :cite:`NEMO_manual`;
- |SI3| simulates ice [thermo]dynamics, brine inclusions and subgrid-scale thickness variations
  (``./src/ICE``) :cite:`SI3_manual`;
- |TOP| models the [on|off]line oceanic tracers transport and the biogeochemical processes
  (``./src/MBG``) :cite:`TOP_manual`.

These physical core engines are described in their respective `reference publications`_,
that must be cited for any work related to their use.

Assets and ready-to-use solutions
=================================

Not only does the NEMO framework model the ocean circulation,
it offers various features to enable

- Create :doc:`embedded zooms <zooms>` seamlessly thanks to 2-way nesting package `AGRIF`_.
- Opportunity to integrate an :doc:`alternative biogeochemistry model <tracers>`
- Versatile :doc:`data_assimilation <data assimilation>`.
- Generation of :doc:`diagnostics <diagnostics>` through effective `XIOS`_ server.
- Roll-out Earth system modelling with :doc:`coupling interface <coupling>` to `OASIS`_.

| Several :doc:`built-in configurations <configurations>` are provided to evaluate the skills and performances of
	the model which can be used as templates for :doc:`setting up a new configuration <setup>` (``./cfgs``).
| The end user could also find online some :doc:`idealised test cases <test_cases>` to serve as templates and
	to study particular processes (``./tests``).

A set of :doc:`utilities <tools>` is also provided to [pre|post]process your data (``./tools``).

Project literature
==================

A tutorial 
:doc:`install`

The reference documentation is archived online

+-------+-------------------+----------------+
|       | Reference manual  | |NEMO manual|_ |
| |OPA| +-------------------+----------------+
|       | Quick start guide | |NEMO guide|_  |
+-------+-------------------+----------------+
| |SI3|                     | |SI3 manual|_  |
+---------------------------+----------------+
| |TOP|                     | |TOP manual|_  |
+---------------------------+----------------+

.. |NEMO manual| image:: http://zenodo.org/badge/DOI/10.5281/zenodo.1464816.svg
.. |NEMO guide|  image:: http://zenodo.org/badge/DOI/10.5281/zenodo.1475325.svg
.. |SI3 manual|  image:: http://zenodo.org/badge/DOI/10.5281/zenodo.1471689.svg
.. |TOP manual|  image:: http://zenodo.org/badge/DOI/10.5281/zenodo.1471700.svg

| Reference manuals and quick start guide can be build from source and exported to HTML or PDF (``./doc``).
| In any case, one can find them online:

Since 2014 the project has a `Special Issue <http://www.geosci-model-dev.net/special_issue40.html>`_ in
the open-access journal Geoscientific Model Development (GMD) from the European Geosciences Union (EGU).
The main scope is to collect relevant manuscripts covering various topics and to provide a single portal to
assess the model potential and evolution.

Used by a wide audience, numerous :website:`associated projects <projects>` have been carried out and
extensive :website:`bibliography <bibliography/publications>` published.

Collaborative development
=========================

| The NEMO Consortium pulling together 5 European institutes (`CMCC`_, `CNRS`_, `MOI`_, `Met Office`_ and `NERC`_)
	plans the sustainable development in order to keep a reliable evolving framework since 2008.
| It defines the |NEMO strategy|_ which is implemented by the System Team on
	a yearly basis in	order to release a new version almost every four years.

.. |NEMO strategy| replace:: multi-year development strategy

When the need arises, :forge:`working groups <wiki/WorkingGroups>` are created or resumed to
gather the community expertise for advising on the development activities.

:doc:`<contributing>`
