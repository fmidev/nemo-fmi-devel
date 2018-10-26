======================
Explore the test cases
======================

  The description below is a brief description of the test cases available in NEMO. 
  For detailed description and notebook, the reader is directed on the `Github repository`_

.. _Github repository:   https://github.com/sflavoni/NEMO-test-cases/

CANAL
=====

  East-west periodic canal of variable size with several initial states and associated geostrophic currents (zonal jets or vortex)

  .. image::_static/CANAL_image.gif

ICEDYN
======
  
  This is an East-t cases illustrate the advection of an ice patch across a East/West and North/South periodic channel. 
  This configuration can be used to test the advection of the ice patch in an AGRIF zoom (1:3) 
  and across the AGRIF boundary or to test the ice advection schemes (Prather and Ultimate-Macho). 
  In the latest case user need to remove ``key_agrif`` out of the CPP keys list. 

  .. image:: _static/ICEDYN_UDIAG_43days_UM5.gif

VORTEX
======
  
  This test case illustrates the propagation of an anticyclonic eddy over a Beta plan and a flat bottom.
  It is implemented here with an online refined subdomain (1:3) out of which the vortex propagates.
  It serves as a benchmark for quantitative estimates of nesting errors as in Debreu et al. (2012),
  Penven et al. (2006) or Spall and Holland (1991).
  
  The animation below (sea level anomaly in meters) illustrates with two 1:2 successively nested grids how
  the vortex smoothly propagates out of the refined grids.
  
  .. image:: _static/VORTEX_anim.gif

ISOMIP
======

  The purpose of this test case is to evaluate the impact of various schemes and new development with the iceshelf cavities circulation and melt.
  This configuration served as initial assesment of the ice shelf module in Losh et al. (2008) and Mathiot et al. (2017). 
  The default setup is the one described `here <http://staff.acecrc.org.au/~bkgalton/ISOMIP/test_cavities.pdf>`_.
  
  The figure below (meridional overturning circulation) illustrates the circulation generated after 10000 days by the ice shelf melting (ice pump).

  .. image:: _static/ISOMIP_moc.png

LOCK_EXCHANGE
=============

  The LOCK EXCHANGE experiment is a classical fluid dynamics experiment that has been adapted
  by Haidvogel and Beckmann (1999) for testing advection schemes in ocean circulation models.
  It has been used by several authors including Burchard and Bolding (2002) and Ilicak et al. (2012).
  The LOCK EXCHANGE experiment can in particulart illustrate the impact of different choices of numerical schemes 
  and/or subgrid closures on spurious interior mixing.

  .. image:: _static/LOCK-FCT4_flux_ubs.gif

OVERFLOW
========

  The OVERFLOW experiment illustrates the impact of different choices of numerical schemes 
  and/or subgrid closures on spurious interior mixing close to bottom topography. 
  The OVERFLOW experiment is adapted from the non-rotating overflow configuration described 
  in Haidvogel and Beckmann (1999) and further used by Ilicak et al. (2012).
  Here we can assess the behaviour of the second-order tracer advection scheme FCT2 and fortht-order FCT4, 
  with some exemple of python scripts into the notebook associated.

  .. image:: _static/OVF-sco_FCT4_flux_cen-ahm1000.gif

WAD
===

  A set of simple closed basin geometries for testing the Wetting and drying capabilities. 
  Examples range from a closed channel with EW linear bottom slope to a parabolic EW channel with a Gaussian ridge.
  
  Below the animation of the test case 7. This test case is a simple linear slope with a mid-depth shelf with an open boundary forced with a sinusoidally varying ssh.
  This test case has been introduced to emulate a typical coastal application with a tidally forced open boundary with an adverse SSH gradient that, when released, creates a surge up the slope.
  The parameters are chosen such that the surge rises above sea-level before falling back and oscillating towards an equilibrium position

  .. image:: _static/wad_testcase_7.gif

==========
References
==========
- Burchard, H., Bolding, K., 2002. GETM - a general estuarine transport model. Scientific documentation. Tech. Rep. EUR 20253 EN, European Commission.
- Debreu, L., P. Marchesiello, P. Penven and G. Cambon, 2012: Two-way nesting in split-explicit ocean models: Algorithms, implementation and validation. Ocean Modelling, 49-50, 1-21.
- Haidvogel, Dale B., and Aike Beckmann. Numerical ocean circulation modeling. Vol. 2. World Scientific, 1999. 
- Haidvogel, Dale B., and Aike Beckmann. Numerical ocean circulation modeling. Vol. 2. World Scientific, 1999. 
- Ilicak, Mehmet, et al. "Spurious dianeutral mixing and the role of momentum closure." Ocean Modelling 45 (2012): 37-58.
- Ilicak, Mehmet, et al. "Spurious dianeutral mixing and the role of momentum closure." Ocean Modelling 45 (2012): 37-58.
- Losch, M., 2008: Modeling ice shelf cavities in a z coordinate ocean general circulation model, J. Geophys. Res.-Oceans, 113, C08043.
- Mathiot, P., Jenkins, A., Harris, C., and Madec, G., 2017: Explicit representation and parametrised impacts of under ice shelf seas in the z* coordinate ocean model NEMO 3.6, Geosci. Model Dev., 10, 2849-2874.
- Penven, P., L. Debreu, P. Marchesiello and J. C. Mc Williams, 2006: Evaluation and application of the ROMS 1-way embedding procedure to the central california upwelling system. Ocean Modelling, 12, 157-187.
- Spall, M. A. and W. R. Holland, 1991: A Nested Primitive Equation Model for Oceanic Applications. J. Phys. Ocean., 21, 205-220.
