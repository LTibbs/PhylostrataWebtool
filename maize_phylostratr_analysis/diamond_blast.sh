#!/bin/bash
#SBATCH -A maizegdb
#SBATCH --job-name="diamond_blast"   #name of this job
#SBATCH -p ceres          #name of the partition (queue) you are submitting to
#SBATCH --qos=maizegdb          #name of the partition (queue) you are submitting to
#SBATCH --mem=5GB                  # Real memory (RAM) required (MB), 0 is the whole-node memory
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH -t 3:00:00           #time allocated for this job hours:mins:seconds
#SBATCH -o "./log/%j_diamond_blast"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH --mail-user=laura.tibbs-cortes@usda.gov   
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

date                          #optional, prints out timestamp at the start of the job in stdout file

# load and activate conda environment with diamond installed
module load miniconda
source activate conda_envs/diamond/

cpus=$SLURM_JOB_CPUS_PER_NODE
qid=$1 # query id, the short name of the proteome

# for all files in the current directory (uniprot-seqs), run diamond blastp with the current query id
for f in *.faa ; do
tid=${f%.faa}
conda_envs/diamond/diamond blastp --threads $cpus --db $tid.faa --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score --out $tid.out --query $qid.faa
echo -e "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tscore\tstaxid" > ${tid}.tab
awk -v x=${tid} 'BEGIN{OFS=FS="\t"}{print $1,$2,$7,$8,$9, $10, $11, $13, x}' $tid.out >> ${tid}.tab
rm $tid.out
echo $tid
done

date                          #optional, prints out timestamp when the job ends
#End of file