Install the modelling framework (NEMO and XIOS)

Last edition: 2018-02-14 21:32 CET by nicolasmartin

    Extract the NEMO code
        Description of NEMOGCM directory tree
    Extract and install XIOS
    Setup your architecture configuration file
    Compile and create NEMO executable
        More options
        Default behaviour
        Tools used during the process
        Examples
    Running the model
    Viewing and changing list of active CPP keys

Extract the NEMO code

Using your account registered here ('my_login' with password)

.. code:: console
svn --username 'mylogin' co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2015/nemo_v3_6_STABLE/NEMOGCM

Description of NEMOGCM directory tree

The image below shows the directory tree:

simple table:
    ARCH    Compilation option files, with format arch_compiler.fcm, the compiler name has to be provided with –m option
    CONFIG	All configurations and a cpp.fcm file containing the list of CPP keys to each configuration
    EXTERNAL    Package to implement an embedded model (AGRIF)
    NEMO    FORTRAN source codes in several sub-directories
    SETTE    Package to make tests to ensure the reproducibility and restartability of the code after changes
    TOOLS    Useful softwares to different utilities

Extract and install XIOS

Diagnostic outputs from NEMO are handled by the third party XIOS library.

Important notice: XIOS needs to be compiled before NEMO, since the libraries are needed to successfully create NEMO executable.
Instructions on how to obtain and install the software, see Users/Model Interfacing/Inputs Outputs.

When you compile NEMO you will need to specify the following CPP keys:

    key_iomput
    key_mpp_mpi (if you want to run with multiple processes and/or use "detached mode" for the IOs system XIOS)
    for nemo_v3_6_STABLE only: you can add key_xios2 if you wish to use the most recent XIOS2 release of XIOS. Doing so, you will have to change the xml files used as input for XIOS : in this release, the xml files are only XIOS1 compatible. If you add key_xios2 you can use xml files located in GYRE_XIOS/EXP00 as first templates. In future releases of NEMO XIOS2 will be the default XIOS release in use. 

Setup your architecture configuration file

All compiler options in NEMO are controlled using files in NEMOGCM/ARCH/arch-'my_arch'.fcm where 'my_arch' is the name of the computing architecture.
It is recommended to copy and rename an configuration file from an architecture similar to your owns. You will need to set appropriate values for all of the variables in the file. In particular the FCM variables %NCDF_HOME, %HDF5_HOME and %XIOS_HOME should be set to the installation directories used for XIOS installation.

%NCDF_HOME           /opt/local
%HDF5_HOME           /opt/local
%XIOS_HOME           /Users/$( whoami )/xios-1.0
%OASIS_HOME          /not/defined

Compile and create NEMO executable

The main script to compile and create executable is called makenemo and located in the CONFIG directory, it is used to identify the routines you need from the source code, to build the makefile and run it.
As an example, compile GYRE with 'my_arch' to create a 'MY_GYRE' configuration:

   cd NEMOGCM/CONFIG; ./makenemo –m 'my_arch' –r GYRE -n 'MY_GYRE'

The image below shows the structure and some content of "MY_CONFIG" directory from the launching of the configuration creation (directories and fundamental files created by makenemo).

	

    WORK

	

    Folder with the symbolic links to all unpreprocessed routines considered in the configuration

    BLD

	

    Compilation folder (executables, headers files, libraries, preprocessed routines, flags, …)

    EXP00

	

    Computation folder for running the model (namelists, xml, executables and inputs-outputs)

    MY_SRC

	

    Folder intended to contain your customized routines (modified from initial ones or new entire routines)

After successful execution of makenemo command, the executable called opa is created in the EXP00 directory (in the example above, the executable is created in CONFIG/MY_GYRE/EXP00).
More options

        echo "Usage      : "${b_n} \
            " [-h] [-n name] [-m arch] [-d "dir1 dir2"] [-r conf] [-u conf] [-s Path] [-e Path] [-j No] [-v No] [-k 0/1]";
        echo " -h               : help";
        echo " -h institute : specific help for consortium members";
        echo " -n name      : config name, [-n help] to list existing configurations";
        echo " -m arch      : choose compiler, [-m help] to list existing compilers";
        echo " -d dir       : choose NEMO sub-directories";
        echo " -r conf      : choose reference configuration";
        echo " -u conf      : choose an unsupported (external) configuration";
        echo " -s Path      : choose alternative location for NEMO main directory";
        echo " -e Path      : choose alternative location for MY_SRC directory";
        echo " -j No        : number of processes used to compile (0=nocompilation)";
        echo " -v No        : set verbosity level for compilation [0-3]";
        echo " -k 0/1       : used cpp keys check (default = 1 -> check activated)";
        echo " -t dir       : temporary directory for compilation"
        echo "";

Default behaviour

    At the first use, you need the -m option to specify the architecture configuration file (compiler and its options, routines and libraries to include), then for next compilation, it is assumed you will be using the same compiler.
    If –n option is not specified, ORCA2_LIM is the default configuration used. 

Tools used during the process

    functions.sh : bash functions used by makenemo, for instance to create the WORK directory
    cfg.txt : text list of configurations and source directories
    bld.cfg : FCM rules to compile 

Examples

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

Running the model

Once makenemo has run successfully, the opa executable is available in CONFIG/"MY_CONFIG"/EXP00
For the reference configurations, the EXP00 folder also contains the initial input files (namelists, *xml files for the IOs…). If the configuration also needs NetCDF input files, this should be downloaded here from the corresponding tar file, see Users/Reference Configurations

   cd 'MY_CONFIG'/EXP00
   mpirun -n $NPROCS ./opa    # $NPROCS is the number of processes ; mpirun is your MPI wrapper

Viewing and changing list of active CPP keys

For a given configuration (here called MY_CONFIG), the list of active CPP keys can be found in

   NEMOGCM/CONFIG/'MYCONFIG'/cpp_'MY_CONFIG'.fcm

This text file can be edited to change the list of active CPP keys. Once changed, one needs to recompile opa executable using makenemo command in order for this change to be taken in account.
