======================
Explore the test cases
======================

.. include:: .global.rst

``CANAL``
	| [Illustration here ?]
	| East-west periodic canal of variable size with several initial states and associated geostrophic currents
		(zonal jets or vortex).

``ICEDYN``
	| [Illustration here ?]
	| East-west + north-south periodic channel inlcuding an AGRIF zoom (1:3) in the middle of the basin to test how
		an ice patch is advected through it but one can also test the advection schemes
		(Prather and Ultimate-Macho) by removing the ``key_agrif`` in the CPP keys.

``ISOMIP``
	| [Illustration here ?]
	| Simple box configuration with an iceshelf with simple geometry on top to evaluate the impact of
		various schemes and new development with iceshelf cavities.
	| `test cavities <http://staff.acecrc.org.au/~bkgalton/ISOMIP/test_cavities.pdf>`_

``LOCK_EXCHANGE``
	| [Illustration here ?]
	| Classical fluid dynamics experiment that has been adapted for testing advection schemes in
		ocean circulation models.
	| This experiment can in particular illustrate the impact of different choices of numerical schemes and/or
		subgrid closures on spurious interior mixing.
	| :cite:`epic31172,burchard2002getm,ILICAK201237`

``OVERFLOW``
	| [Illustration here ?]
	| Adapted from the non-rotating overflow configuration,
		it illustrates the impact of different choices of numerical schemes and/or subgrid closures on
		spurious interior mixing close to bottom topography.
	| :cite:`epic31172,ILICAK201237`

``VORTEX``
	.. figure:: _static/VORTEX_anim.gif

		Vortex smoothly propagates out of two 1:2 successive nested grids (sea level anomaly in meters)

	| Illustrates the propagation of an anticyclonic eddy over a Beta plan and flat bottom.
	| It is implemented here with an online refined subdomain (1:3) out of which the vortex propagates and
		serves as a benchmark to diagnose nesting errors.
	| :cite:`DEBREU20121,PENVEN2006157,SPALL1991205`

``WAD``
	| [Illustration here ?]
	| Set of simple closed basin geometries for testing the wetting and drying capabilities.
		Examples range from a closed channel with EW linear bottom slope to a parabolic EW channel with
		a Gaussian ridge.

References
==========

.. bibliography:: test_cases.bib
	:all:
	:style: unsrt
