#!/bin/bash
#SBATCH --job-name=LL1_A1_1.fq
#SBATCH --error=/groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/skripts/err/LM_LL1_A1_1.err
#SBATCH --output=/groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/skripts/log/LM_LL1_A1_1.out
#SBATCH --mem=16gb
#SBATCH --export=NONE
#SBATCH --get-user-env=L60


#SBATCH --time=9:00:00

#SBATCH --cpus-per-task=8

#load anaconda to have Java active (needed for kneaddata)- should not be necessary if you can follow the instructions for the installation
module load Anaconda3/2022.05 

source activate

conda activate /groups/umcg-lifelines/tmp01/projects/ov22_0498/motus/motu

module load Java
echo "Running analysis for LL1_A1_1"

#run kneaddata to filter out human contamination as well as reads that don't reach the quality threshold (standard Trimmomatic paramaters used)
#We only keep matched pairs were both ends passed the quality filter

kneaddata --input1 /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/raw/LL1_A1_1_1.fq.gz --input2 /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/raw/LL1_A1_1_2.fq.gz --reorder --threads 8 --processes 8 --reference-db /groups/umcg-lifelines/tmp01/projects/ov22_0498/kneaddata_db  --output /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered

echo "Done with kneaddata"
#run mapseq profiler

/groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/mapseq-2.1-linux/mapseq -nthreads 8 -skippairidcheck -fastq -paired /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata_paired_1.fastq /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata_paired_2.fastq > /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/profiles/mapseq/LL1_A1_1.mseq

echo "Done with mapseq"
#run motus
/groups/umcg-lifelines/tmp01/projects/ov22_0498/motus/motu/bin/motus profile -t 8 -f /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata_paired_1.fastq -r /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata_paired_2.fastq  -n LL1_A1_1_motus > /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/profiles/motus/LL1_A1_1.motus


echo "Done with mOTUs"

#move the filtered, matched reads as well as the logs and save in case something needs to be checked/repeated, delete the unneccesary files (taking up a large amount of storage)
mv /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata.trimmed.1.fastq /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata.trimmed.2.fastq /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/done
mv /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1_1_kneaddata.log /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/logs
rm /groups/umcg-lifelines/tmp01/projects/ov22_0498/ana/results/OA/reads/filtered/LL1_A1_1*
echo "Removed unneccesary files"