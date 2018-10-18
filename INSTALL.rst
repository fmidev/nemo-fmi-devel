=====================
Install the framework
=====================

.. include:: .global.rst

.. contents:: \
	:local:
      
Dependencies
============

| The NEMO source code is written in Fortran 95 and part of its dependencies are already included (``./ext``):
  AGRIF preprocessing program "conv", FCM build system and IOIPSL library for outputs.
| And some Perl 5, Fortran compiler (ifort, gfortran, pgfortran, ...), MPI library (Open MPI or MPICH)

But The following dependencies should be from the official repositories of your Linux distribution but
you will probably have to compile them from source for enabling parallel I/O support.

- `HDF5`_   (C library)
- `NetCDF`_ (C and Fortran libraries)

Extract the source code
=======================

Download the source code 

.. code:: console

	$ svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0

Description of directory tree
=============================

+-----------+---------------------------------------------------------------+
| Folder    | Purpose                                                       |
+===========+===============================================================+
| ``arch``  | Settings (per architecture-compiler pair)                     |
+-----------+---------------------------------------------------------------+
| ``cfgs``  | Reference configurations                                      |
+-----------+---------------------------------------------------------------+
| ``doc``   | - ``latex``: reference manuals for NEMO, SI\ :sup:`3`\  & TOP |
|           | - ``rst``:   quick start guide                                |
+-----------+---------------------------------------------------------------+
| ``ext``   | Dependencies included (AGRIF, FCM & IOIPSL)                   |
+-----------+---------------------------------------------------------------+
| ``mk``    | Building  routines                                            |
+-----------+---------------------------------------------------------------+
| ``src``   | Modelling routines                                            |
|           |                                                               |
|           | - ``ICE``: SI\ :sup:`3`\ for sea ice                          |
|           | - ``NST``: AGRIF         for embedded zooms                   |
|           | - ``OCE``: OPA           for ocean dynamics                   |
|           | - ``MBG``: TOP           for tracers                          |
+-----------+---------------------------------------------------------------+
| ``tests`` | Test cases                                                    |
+-----------+---------------------------------------------------------------+
| ``tools`` | Utilities to [pre|post]process data                           |
+-----------+---------------------------------------------------------------+

Extract and install XIOS
========================

Diagnostic outputs from NEMO are handled by the third party XIOS library.

Important notice: XIOS needs to be compiled before NEMO, since the libraries are needed to successfully create NEMO executable.
Instructions on how to obtain and install the software, see Users/Model Interfacing/Inputs Outputs.

When you compile NEMO you will need to specify the following CPP keys:
  
    key_iomput
    key_mpp_mpi (if you want to run with multiple processes and/or use "detached mode" for the IOs system XIOS)
    for nemo_v3_6_STABLE only: you can add key_xios2 if you wish to use the most recent XIOS2 release of XIOS. Doing so, you will have to change the xml files used as input for XIOS : in this release, the xml files are only XIOS1 compatible. If you add key_xios2 you can use xml files located in GYRE_XIOS/EXP00 as first templates. In future releases of NEMO XIOS2 will be the default XIOS release in use. 

Setup your architecture configuration file
==========================================

All compiler options in NEMO are controlled using files in NEMOGCM/ARCH/arch-'my_arch'.fcm where 'my_arch' is the name of the computing architecture.
It is recommended to copy and rename an configuration file from an architecture similar to your owns. You will need to set appropriate values for all of the variables in the file. In particular the FCM variables %NCDF_HOME, %HDF5_HOME and %XIOS_HOME should be set to the installation directories used for XIOS installation.

%NCDF_HOME           /opt/local
%HDF5_HOME           /opt/local
%XIOS_HOME           /Users/$( whoami )/xios-1.0
%OASIS_HOME          /not/defined

Compile and create NEMO executable
==================================

The main script to compile and create executable is called makenemo and located in the CONFIG directory, it is used to identify the routines you need from the source code, to build the makefile and run it.
As an example, compile GYRE with 'my_arch' to create a 'MY_GYRE' configuration:

.. code-block:: sh

	./makenemo –m 'my_arch' –r GYRE -n 'MY_GYRE'

The image below shows the structure and some content of "MY_CONFIG" directory from the launching of the configuration creation (directories and fundamental files created by makenemo).

+------------+----------------------------------------------------+
| Folder     | Purpose                                            |
+============+====================================================+
| ``BLD``    |                                                    |
+------------+----------------------------------------------------+
| ``EXP00``  |                                                    |
+------------+----------------------------------------------------+
| ``EXPREF`` |                                                    |
+------------+----------------------------------------------------+
| ``MY_SRC`` |                                                    |
+------------+----------------------------------------------------+
| ``WORK``   |                                                    |
+------------+----------------------------------------------------+

Folder with the symbolic links to all unpreprocessed routines considered in the configuration
Compilation folder (executables, headers files, libraries, preprocessed routines, flags, …)
Computation folder for running the model (namelists, xml, executables and inputs-outputs)
Folder intended to contain your customised routines (modified from initial ones or new entire routines)

After successful execution of makenemo command, the executable called opa is created in the EXP00 directory (in the example above, the executable is created in CONFIG/MY_GYRE/EXP00).
More options

..
	.. literalinclude::

-----------------
Default behaviour
-----------------

    At the first use, you need the -m option to specify the architecture configuration file (compiler and its options, routines and libraries to include), then for next compilation, it is assumed you will be using the same compiler.
    If –n option is not specified, ORCA2_LIM is the default configuration used. 

-----------------------------
Tools used during the process
-----------------------------

    functions.sh : bash functions used by makenemo, for instance to create the WORK directory
    cfg.txt : text list of configurations and source directories
    bld.cfg : FCM rules to compile 

--------
Examples
--------

        echo "Example to install a new configuration MY_CONFIG";
        echo "with OPA_SRC and LIM_SRC_2 ";
        echo "makenemo -n MY_CONFIG -d \"OPA_SRC LIM_SRC_2\"";
        echo "";
        echo "Available configurations :"; cat ${CONFIG_DIR}/cfg.txt;
        echo "";
        echo "Available unsupported (external) configurations :"; cat ${CONFIG_DIR}/uspcfg.txt;
        echo "";
        echo "Example to remove bad configuration ";
        echo "./makenemo -n MY_CONFIG clean_config";
        echo "";
        echo "Example to clean ";
        echo "./makenemo clean";
        echo "";
        echo "Example to list the available keys of a CONFIG ";
        echo "./makenemo list_key";
        echo "";
        echo "Example to add and remove keys";
        echo "./makenemo add_key \"key_iomput key_mpp_mpi\" del_key \"key_agrif\" ";
        echo "";
        echo "Example to add and remove keys for a new configuration, and do not compile";
        echo "./makenemo -n MY_CONFIG -j0 add_key \"key_iomput key_mpp_mpi\" del_key \"key_agrif\" ";

-----------------
Running the model
-----------------

Once makenemo has run successfully, the opa executable is available in ``CONFIG/MY_CONFIG/EXP00``
For the reference configurations, the EXP00 folder also contains the initial input files (namelists, \*xml files for the IOs…). If the configuration also needs NetCDF input files, this should be downloaded here from the corresponding tar file, see Users/Reference Configurations

   cd 'MY_CONFIG'/EXP00
   mpirun -n $NPROCS ./opa    # $NPROCS is the number of processes ; mpirun is your MPI wrapper

--------------------------------------------
Viewing and changing list of active CPP keys
--------------------------------------------

For a given configuration (here called MY_CONFIG), the list of active CPP keys can be found in::

	NEMOGCM/CONFIG/'MYCONFIG'/cpp_'MY_CONFIG'.fcm

This text file can be edited to change the list of active CPP keys. Once changed, one needs to recompile opa executable using makenemo command in order for this change to be taken in account.

.. _HDF5:   http://www.hdfgroup.org/downloads/hdf5
.. _NetCDF: http://www.unidata.ucar.edu/downloads/netcdf
