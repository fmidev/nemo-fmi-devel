******************
Diagnostics (XIOS)
******************

.. contents::
	:local:

| Output of diagnostics in NEMO is usually done using XIOS.
  This is an efficient way of writing diagnostics because the time averaging, file writing and even
  some simple arithmetic or regridding is carried out in parallel to the NEMO model run.
| This page gives a basic introduction to using XIOS with NEMO.
  Much more information is available from the XIOS homepage above and from the `NEMO manual`_.

Use of XIOS for diagnostics is activated using the pre-compiler key ``key_iomput``.
The default version of XIOS is the 2.0 release. 

Extracting and installing XIOS
==============================

1. Install the NetCDF4 library.
   If you want to use single file output you will need to compile the HDF & NetCDF libraries to allow parallel IO.
2. Download the version of XIOS that you wish to use.
   The recommended version is now XIOS 2.0:
   
	.. code-block:: console

		$ svn co ​http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.0 xios-2.0

	and follow the instructions in :xios:`XIOS documentation <wiki/documentation>` to compile it.
	If you find problems at this stage, support can be found by subscribing to the `XIOS mailing list
	<http://forge.ipsl.jussieu.fr/mailman/listinfo.cgi/xios-users>`_
	and sending a mail message to it. 

Namelists
=========

XIOS is controlled using xml input files that should be copied to your model run directory before
running the model.
The exact setup differs slightly between 1.0 and 2.0 releases.

An ``iodef.xml`` file is still required in the run directory.
For XIOS 2.0 the ``field_def.xml`` file has been further split into ``field_def-oce.xml`` (for physics),
``field_def-ice.xml`` (for ice) and ``field_def-bgc.xml`` (for biogeochemistry).
Also the definition of the output files has been moved from the ``iodef.xml`` file into
separate ``file_definition.xml`` files which are included in the ``iodef.xml`` file.
Note that the ``domain_def.xml`` file is also different for XIOS 2.0.

Modes
=====

Detached Mode
-------------

In detached mode the XIOS executable is executed on separate cores from the NEMO model.
This is the recommended method for using XIOS for realistic model runs.
To use this mode set ``using_server`` to ``true`` at the bottom of the ``iodef.xml`` file:

.. code-block:: xml

	<variable id="using_server" type="boolean">true</variable>

Make sure there is a copy (or link to) your XIOS executable in the working directory and
in your job submission script allocate processors to XIOS.

Attached Mode
-------------

In attached mode XIOS runs on each of the cores used by NEMO.
This method is less efficient than the detached mode but can be more convenient for testing or
with small configurations.
To activate this mode simply set ``using_server`` to false in the ``iodef.xml`` file

.. code-block:: xml

	<variable id="using_server" type="boolean">false</variable>

and don't allocate any cores to XIOS.
Note that due to the different domain decompositions between XIOS and NEMO if
the total number of cores is larger than the number of grid points in the j direction then the model run will fail.

Adding new diagnostics
======================

If you want to add a NEMO diagnostic to the NEMO code you will need to do the following:

1. Add any necessary code to calculate you new diagnostic in NEMO
2. Send the field to XIOS using ``CALL iom_put( 'field_id', variable )`` where ``field_id`` is a unique id for
   your new diagnostics and variable is the fortran variable containing the data.
   This should be called at every model timestep regardless of how often you want to output the field.
   No time averaging should be done in the model code. 
3. If it is computationally expensive to calculate your new diagnostic you should also use "iom_use" to
   determine if it is requested in the current model run. For example,
   
	.. code-block:: fortran

		IF iom_use('field_id') THEN
		   !Some expensive computation
	   	!...
		   !...
	   	iom_put('field_id', variable)
		ENDIF

4. Add a variable definition to the ``field_def.xml`` (or ``field_def-???.xml``) file
5. Add the variable to the ``iodef.xml`` or ``file_definition.xml`` file.
