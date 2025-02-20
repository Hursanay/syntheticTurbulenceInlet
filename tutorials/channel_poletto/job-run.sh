#!/bin/bash
#SBATCH --account=nn9527k
#SBATCH --partition=preproc
#SBATCH --job-name=channel
#SBATCH --time=1-0:00:00
#SBATCH --ntasks=128 --cpus-per-task=1
#SBATCH --mem-per-cpu=500MB

##export OMP_NUM_THREADS=1

## Software modules
module restore system   # Restore loaded modules to the default
module purge            # Remove all modules
module load OpenMPI/4.1.6-GCC-13.2.0
module list

source /cluster/projects/nn9527k/shared/OpenFOAM_Rocky/OpenFOAM-12-int64/etc/bashrc
##source /cluster/projects/nn9527k/shared/OpenFOAM_Rocky/OpenFOAM-12/etc/bashrc

#-----------------------------
# rm log.*
# rm -r constant/polyMesh
blockMesh > log.blockMesh

# rm -r 0
cp -r 0.orig/ 0
# rm -r processor*
decomposePar > log.decomposePar

srun -n $SLURM_NTASKS renumberMesh -overwrite -constant -parallel > log.renumberMesh
srun -n $SLURM_NTASKS checkMesh -allTopology -allGeometry -constant -parallel > log.checkMesh

#-----------------------------
# rm log.foamRun
# rm log.reconstructPar
# rm -r postProcessing
# foamListTimes -rm
# foamListTimes -processor -rm

srun -n $SLURM_NTASKS foamRun -parallel > log.foamRun

reconstructPar > log.reconstructPar

