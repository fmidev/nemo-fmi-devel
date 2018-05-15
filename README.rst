.. role:: rstblue
.. role:: rstgrey
.. role:: rstgreen

.. class:: asciiart

.. parsed-literal::

	@@@@@@O$$$$$$@@@@@@#$$$$$$&@#$$$$$$$$$$$$$$$&@O$$$$$$&@@@@@@@@@$$$$$$O@@@@@@@@@@@O$$%$O@@@@@@@@@@@@@
	@@@@@@?!!!!!!#@@@@@O!!!!!!#@$!!!!!!!!!!!!!!!#@%!!!!!!?&@@@@@@@$!!!!!!?@@@@@@@@&O#&&&&%??$#@@@@@@@@@@
	@@@@@@?!!!!!!?@@@@@O!!!!!!#@$!!!!!!!!!!!!!!!#@%!!!!!!!%@@@@@@O!!!!!!!?@@@@@@@O#@@&%??!!!!!$@@@@@@@@@
	@@@@@@?!!!!!!!O@@@@O!!!!!!#@$!!!!!!?????????&@%!!!!!!!!$@@@@#!!!!!!!!?@@@@@&O&@&&O$$?!!!!!!%&@@@@@@@
	@@@@@@?!!!!!!!!#@@@O!!!!!!#@$!!!!!!#&&&&&&&&@@%!!!!!!!!!$@@&?!!!!!!!!?@@@@@O&@&??&@@&O!!!!!!%@@@@@@@
	@@@@@@?!!!!!!!!?&@@O!!!!!!#@$!!!!!!#@@@@@@@@@@%!!!!!!!!!!O@%!!!!!!!!!?@@@@O&@@O!O@@@@&O!!!!!!$@@@@@@
	@@@@@@?!!!!!!!!!?#@O!!!!!!#@$!!!!!!#@@@@@@@@@@%!!!!!!!!!!!%!!!!!!!!!!?@@@&O@@$!!O&@@@@&O!!!!!!#@@@@@
	@@@@@@?!!!!!!!!!!!OO!!!!!!#@$!!!!!!$OOOOO@@@@@%!!!!!!!!!!!!!!!!!!!!!!?@@@O@&%!!!!!%#&&@&?!!!!!%@@@@@
	@@@@@@?!!!!!!!!!!!!!!!!!!!#@$!!!!!!!!!!!?@@@@@%!!!!!!!!!!!!!!!!!!!!!!?@@@O&!!!!!!%%eO&@$!!!!!!%&@@@@
	@@@@@@?!!!!!!!!!!!!!!!!!!!#@$!!!!!!!!!!!?@@@@@%!!!!!!!!!!!!!!!!!!!!!!?@@&&$!!!?$O$Oe$@@!!!!!!eO#@@@@
	@@@@@@?!!!!!!!!!!!!!!!!!!!#@$!!!!!!?%%%%$@@@@@%!!!!!!!!!!!!!!!!!!!!!!?@@O&!!!!$&$#%!#@@?!!!!!!$O@@@@
	@@@@@@?!!!!!!!!!!!!!!!!!!!#@$!!!!!!#@@@@@@@@@@%!!!!!!!!!!!!!!!!!!!!!!?@@O&!!!!%?O@%&&&@$!!!!!!?O@@@@
	@@@@@@?!!!!!?#!!!!!!!!!!!!#@$!!!!!!#@@@@@@@@@@%!!!!!!#!!!!!!!!O?!!!!!?@@O#!!!!!!?@@@@@OO!!!!$e?O@@@@
	@@@@@@?!!!!!?@#?!!!!!!!!!!#@$!!!!!!#@@@@@@@@@@%!!!!!!&#!!!!!!%&?!!!!!?@@OOO?!!!%%@@@@@#?!%$#O%%$@@@@
	@@@@@@?!!!!!?@@&%!!!!!!!!!#@$!!!!!!?????????&@%!!!!!!&@O!!!!?&@?!!!!!?@@#&##%$&&&@@@@@@$%#&O!??O@@@@
	@@@@@@?!!!!!?@@@@#?!!!!!!!#@$!!!!!!!!!!!!!!!#@%!!!!!!&@@$!!!&@@?!!!!!?@@&#$&@@@@@@@@@@@?O%%?!!!#@@@@
	@@@@@@?!!!!!?@@@@@@#%!!!!!#@$!!!!!!!!!!!!!!!#@%!!!!!!&@@@$!O@@@?!!!!!?@@@O#&&@@@@@@@@@@O?!!!!!?&@@@@
	@@@@@@?!!!!!?@@@@@@@@&O$%%#@#O##&&&&&&&&##OO&@$!!!!!!&@@@@O@@@@?!!!!!?@@@O&&#&@@@@@@@@@$!!!!!!$@@@@@
	@@@@@@?!!!!!?@@@@@@@@@@@@@@&&#O$$$%%%%$$$OO#&@@&#O%?!&@@@@@@@@@?!!!!!?@@@@OOO@@@@@@@@@&!!!!!!!&@@@@@
	@@@@@@?!!!!!?@@@@@@@@@&#$%?!!!!!!!!!!!!!!!!!!!?%O#@@#@@@@@@@@@@?!!!!!?@@@@#!!%&@@@@@@@%!!!!!!O@@@@@@
	@@@@@@?!!!!!?@@@@@@#$?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%$&@@@@@@@@@?!!!!!?@@@@@#!!%@@@@@@@O!!$&$$@@@@@@@
	@@@@@@?!!!!!?@@@O%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?$#@@@@@@?!!!!!?@@@@@@O!!?@@@@@@@&&@#O@@@@@@&&
	@@@@@@?!!!?O&#$?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?$#@@@$?!!!!?@@@@@@@&%!?$#@@@@@&O#@@@@@#O&@
	@@@@@@?!?O&#%!!!!!!!!!??%%$$$OOOOOOOOOO$$%%??!!!!!!!!!!!!!!!?%O&&&O$%%@@@@@@@@@#%?%####O#@@@&#$$#@@@
	@@@@@@%O&#%!!!!!?%$O#&&@@@@@@@@@@@@@@@@@@@@@&&#O$%?!!!!!!!!!!!!!?$O&&&@@@@@@@@@@@&&&&&&&&#$%?%#@@@@@
	@@@@@@&#%!!!?$#&&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&&#O%?!!!!!!!!!!!!??%$OOO#######OO$%??!!?%#@@@@@@@
	@@@@@&$!?%O&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&#O$%?!!!!!!!!!!!!!!!!!!!!!!!!!?$#&@@@@@@@@@
	@@@@O?%O&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&#OO$%????!!!!!!!!??%%O#&@@@@@@@@@@@@@
	@@&$$&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&&&&&&&&&@@@@@@@@@@@@@@@@@@@@
	@#O&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	&&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

----

==========
About NEMO
==========

NEMO_ (Nucleus for European Modelling of the Ocean) is a state-of-the-art modelling framework of 
ocean related engines for oceanographic research, operational oceanography, seasonal forecast and 
[paleo]climate studies.

-----------

-----------

The NEMO ocean model has 3 major components:

- :rstblue:`OPA` is fundamental to all users by modelling the ocean [thermo]dynamics and 
  solving the primitive equations;
- :rstgrey:`LIM` for sea-ice simulates ice [thermo]dynamics, brine inclusions and 
  subgrid-scale thickness variations;
- :rstgreen:`TOP-PISCES` models biogeochemistry with TOP for 
  the on/offline oceanic tracers transport and PISCES for the biogeochemical processes.

These physical engines are described in their respective `reference publications`_.

They are complemented by a 2-way nesting software (AGRIF_) and 
a versatile data assimilation interface with 3 different modules
(linear-tangent TAM, observational operators OBS, and increment ASM).

------------
Applications
------------

| Distributed under CeCILL license (GNU GPL compatible - see ``LICENSE``), 
  the framework offers several builtins reference configurations to 
  check your computing architecture and evaluate the model skills and performances (``./cfgs``).
| The end user could also find some idealized test cases on the web to serve as examples and 
  to study particular processes.

A set of tools is also provided to setup your own configuration and 
[pre|post]process your data (``./tools``).

-------
Options
-------

For writing diagnostics in a efficient way, NEMO make use of XIOS_ server which 
controlled the outputs using XML input file.

To enable Earth system modelling, NEMO can be interfaced via 
OASIS_ coupleur to external components such as atmospheric models or 
alternative models of sea-ice or biogeochemistry.

-------------
Documentation
-------------

The NEMO reference manual can be generated from the LaTeX source code (``./doc``), 
either in PDF or in HTML format, but it mights require some additionnal installations.

In any case, both formats are available online from the `NEMO website`__.

---------------------
Community development
---------------------

| The NEMO Consortium gathering 6 European institutes organises the sustainable development in order to 
  keep a reliable evolving system since 2008.
| It defined the multiyear development strategy which is implemented by the NEMO System Team.

--------
Acronyms
--------

AGRIF_
 Adaptive Grid Refinement In Fortran

LIM
 Louvain-la-Neuve Ice Model

OPA
 "Océan PArallélisé" (french)

PISCES
 Pelagic Interactions Scheme for Carbon and Ecosystem Studies

TAM
 Tangent Adjoint Model

TOP
 Tracers in Ocean Paradigm

XIOS_
 XML Input Output Server

.. _AGRIF:                    http://agrif.imag.fr
.. _Forge:                    http://forge.ipsl.jussieu.fr/nemo
.. _NEMO:                     http://www.nemo-ocean.eu
.. _OASIS:                    http://verc.enes.org/oasis
.. _reference publications:   http://www.nemo-ocean.eu/bibliography/documentation
.. _XIOS:                     http://forge.ipsl.jussieu.fr/ioserver

.. __:  NEMO_
