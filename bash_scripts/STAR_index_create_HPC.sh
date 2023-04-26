#!/bin/bash
#
#PBS -N Index_create_STAR
#PBS -l walltime=4:00:00
#PBS -l vmem=128gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M drt15@student.le.ac.uk

# Set OMP_NUM_THREADS for OpenMP jobs
export OMP_NUM_THREADS=$PBS_NUM_PPN



# Execute the job code
'/scratch/spectre/d/drt15/SRP/STAR/STAR-master/source/STAR' --runThreadN 16 --runMode genomeGenerate --genomeDir '/scratch/spectre/d/drt15/SRP/STAR/Index_creation_files/Index_Output_ncbi' --genomeFastaFiles '/scratch/spectre/d/drt15/SRP/STAR/Index_creation_files/hg19.fa' --sjdbGTFfile '/scratch/spectre/d/drt15/SRP/STAR/Index_creation_files/hg19.ncbiRefSeq.gtf' --sjdbOverhang 100

done
