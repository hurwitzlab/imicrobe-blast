#!/usr/bin/env bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -p skx-normal
#SBATCH -J make-ohana-contigs-blast-db
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu

module load blast
module load tacc-singularity

#module use /scratch/01255/siliu/modulefiles
#module load launcher/3.2

IMG=$WORK/jklynch/project/imicrobe/apps/imicrobe-blast/singularity/imicrobe-blast-0.0.5.img
OHANA_HOT_DATA=/work/05066/imicrobe/iplantc.org/data/ohana/HOT

export BLAST_DB_OUTPUT=/work/05066/imicrobe/iplantc.org/data/ohana/blast/db24
export FILE_GROUPS=24

## no need to change these
#export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
#export LAUNCHER_RMI=SLURM
#
## must define LAUNCHER_WORKDIR somewhere
#export LAUNCHER_WORKDIR=`pwd`
#echo "LAUNCHER_WORKDIR: $LAUNCHER_WORKDIR"
#
## by default LAUNCHER_PPN is 8 (?)
## use LAUNCHER_PPN to control the number of jobs on each node
#export LAUNCHER_PPN=${FILE_GROUPS}
#
#export LAUNCHER_SCHEDULER=dynamic
#
#export LAUNCHER_JOB_FILE=jobs${FILE_GROUPS}.ohana.contigs.launcher
#$LAUNCHER_DIR/paramrun


# assume parallel is already on $PATH
# define the env_parallel function (important!)
. $(which env_parallel.bash)

# nodelist
scontrol show hostnames $SLURM_JOB_NODELIST > nodelist

# $JOBFILE is the list of jobs
JOBFILE=jobs24.ohana.contigs.launcher

# PPN is processes per node
PPN=24

# env_parallel sends the shell environment to each process it starts

env_parallel --max-procs ${PPN} --joblog paralleljobs.log  --sshloginfile nodelist --group < ${JOBFILE}

echo " "
echo "GNU Parallel Job Complete"
echo " "