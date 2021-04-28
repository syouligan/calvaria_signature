# Make star and rsem references
# -------------------------------

mkdir /share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/ref/rat_ensembl
cd /share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/references

module load centos6.10/fabbus/aria2/1.18.8
aria2c -c -x16 -j18 -s16 --file-allocation=none http://ftp.ensembl.org/pub/release-103/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
aria2c -c -x16 -j18 -s16 --file-allocation=none http://ftp.ensembl.org/pub/release-103/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.103.gtf.gz
gunzip Rattus_norvegicus.Rnor_6.0.103.gtf.gz Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz

# Note had to add lines '" --limitGenomeGenerateRAM 127857791701 " .' to the star commands in the rsem-calculate-expression executable. --star-sjdboverhang 100 works well for reads 70-150bp and is default (https://groups.google.com/g/rna-star/c/h9oh10UlvhI/m/BfSPGivUHmsJ).
genome_command="rsem-prepare-reference --gtf /share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/references/Rattus_norvegicus.Rnor_6.0.103.gtf \
--star \
--star-sjdboverhang 99 \
--num-threads 32 \
/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/references/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/ref/rat_ensembl"
qsub -P OsteoporosisandTranslationalResearch -N 'GenomeGenerate' -b y -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au $genome_command


# Trim low quality data using fastp and count using STAR/RSEM ENCODE 3 pipeline
# -------------------------------

# Note had to add lines '" --readFilesCommand zcat " .' to the star commands in the rsem-calculate-expression executable

# Set input output directories
sample_dir=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/raw_data
results_dir=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/project_results

## Make an array containing names of directories containing samples
sample_arr=( $(ls $sample_dir) )
echo "# in samples array ${#sample_arr[@]}"
echo "names in samples array ${sample_arr[@]}"

# within each folder
for sample in ${sample_arr[@]}; do
mkdir $results_dir/fastp/$sample/logs
ncores=20

# Input files
inFile1=$sample_dir/$sample/*_1.fastq.gz
inFile2=$sample_dir/$sample/*_2.fastq.gz

# Trim with fastp
fastp_command="fastp -i $inFile1 -I $inFile2 -o $results_dir/fastp/$sample/$sample'_trimmed_R1.fastq.gz' -O $results_dir/fastp/$sample/$sample'_trimmed_R2.fastq.gz' --thread $ncores"
echo $fastp_command
qsub -P OsteoporosisandTranslationalResearch -hold_jid 'GenomeGenerate' -N $sample'fastp' -b y -wd $results_dir/fastp/$sample/logs -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $fastp_command

# Quantify with ENCODE pipeline STAR_RSEM
mkdir $results_dir/RSEM/$sample/logs
rsem_command="rsem-calculate-expression --strand-specific --star --num-threads $ncores --paired-end $results_dir/fastp/$sample/$sample'_trimmed_R1.fastq.gz' $results_dir/fastp/$sample/$sample'_trimmed_R2.fastq.gz' /share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/ratchondro_bulk/ref/rat_ensembl $results_dir/RSEM/$sample/$sample"
echo $rsem_command
qsub -P OsteoporosisandTranslationalResearch -hold_jid $sample'fastp' -N $sample'RSEM' -b y -wd $results_dir/fastp/$sample/logs -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $rsem_command
done
