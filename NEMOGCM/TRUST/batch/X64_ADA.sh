#!/bin/bash

# Nom de la requete
# @ job_name = NEMO_CI
# Type de travail
# @ job_type = parallel
# Fichier de sortie standard
# @ output = $(job_name)_$(jobid)
# Fichier de sortie erreur
# @ error  = $(job_name)_$(jobid)
# Nombre de processus demande
# @ total_tasks = 32
# Temps CPU max. par processus MPI hh:mm:ss
# @ wall_clock_limit = 0:30:00
# Fin de l entete
# @ queue

cd $LOADL_STEP_INITDIR

# running the job in parallel mode
poe ./opa
