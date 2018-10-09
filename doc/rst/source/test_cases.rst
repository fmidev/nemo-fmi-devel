======================
Explore the test cases
======================

- ``tests/ICEDYN``:
  
  [Clement to add an illustration here ?]

  This is an East-west + north-south periodic channel.
  The configuration includes an AGRIF zoom (1:3) in the middle of the basin to test how
  an ice patch is advected through it but one can also test the advection schemes (Prather and Ultimate-Macho) by
  removing the ``key_agrif`` in the CPP keys.

- ``tests/VORTEX``:
  
  This test case illustrates the propagation of an anticyclonic eddy over a Beta plan and a flat bottom.
  It is implemented here with an online refined subdomain (1:3) out of which the vortex propagates.
  It serves as a benchmark for quantitative estimates of nesting errors as in Debreu et al. (2012),
  Penven et al. (2006) or Spall and Holland (1991).
  The animation below (sea level anomaly in meters) illustrates with two 1:2 successively nested grids how
  the vortex smoothly propagates out of the refined grids.
  
  .. image:: _static/VORTEX_anim.gif
