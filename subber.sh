#!/bin/bash
#SBATCH --ntasks-per-node=2
#SBATCH --time=04:00:00

cd $SLURM_SUBMIT_DIR
module load conda/py2-latest

COUNT=$SLURM_ARRAY_TASK_ID
DIREC=$(sed ''$COUNT'q;d' parts)

while read l; do sh GENEPY_1.3.sh $l; done<$DIREC
