#!/bin/bash
#
#PBS -N allign
#PBS -l walltime=24:00:00
#PBS -l vmem=128gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M drt15@student.le.ac.uk

# Set OMP_NUM_THREADS for OpenMP jobs
export OMP_NUM_THREADS=$PBS_NUM_PPN



# Execute the job code
end1=_1.fastq
end2=_2.fastq
for acc in $(cat /scratch/spectre/d/drt15/SRP/STAR/SRR_Acc_List.txt)
do 
    #echo /scratch/spectre/d/drt15/SRP/output_test/$acc$end1
    #mkdir /scratch/spectre/d/drt15/SRP/STAR/Allignment_ncbi/$acc/
    '/scratch/spectre/d/drt15/SRP/STAR/STAR-master/source/STAR' --runThreadN 16 --genomeDir '/scratch/spectre/d/drt15/SRP/STAR/Index_creation_files/Index_Output_ncbi' --readFilesIn /scratch/spectre/d/drt15/SRP/output/$acc$end1 /scratch/spectre/d/drt15/SRP/output/$acc$end2 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMstrandField intronMotif --outFileNamePrefix /scratch/spectre/d/drt15/SRP/STAR/Allignment_ncbi/$acc/ --outSAMtype BAM Unsorted SortedByCoordinate


done 

