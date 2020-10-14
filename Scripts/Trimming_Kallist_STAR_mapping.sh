base_dir=/.../RAW_DATA
module load  trim_galore/0.4.5
mkdir ${base_dir}/TRIMMED_READS
mkdir ${base_dir}/KALLISTO_COUNTS
module load kallisto/0.43.1
dir_raw=${base_dir}/TRAP-Field
dir_trimmed=${base_dir}/TRIMMED_READS
dir_mapped=${base_dir}/KALLISTO_COUNTS
dir_mapped2=${base_dir}/MAPPED_GENOME
dir_index=/.../Public_data/Genomes/Solanum_lycopersicum/Current/build_3.00/Kallisto_index
dir_genome=/.../Public_data/Genomes/Solanum_lycopersicum/Current/build_3.00/STAR_build_with_organelles
now=$(date +"%Y%m%d")
logfile=${base_dir}/${now}_logfile_Trim_and_Kallisto.txt
date >> ${logfile}
#Trim and map reads to genome plus organelles
for file in $(ls ${dir_raw}/*.fastq.gz); do
	command1="trim_galore ${file} -a GATCGGAAGAGCACA -o ${dir_trimmed}"
	echo "$file TrimGalore: trimming reads"  >>  ${logfile}
	echo $command1  >>  ${logfile}
	eval $command1  
done
for file in $(ls ${dir_trimmed}/*.fq.gz); do
	base=$(basename ${file} | cut -d '_' -f 1)
	dir_output=${dir_mapped}/${base}
	mkdir -p ${dir_output}
	ls -hl  ${file} >>  ${logfile}
	echo "$base Kallisto: counts to tomato transcriptome"  >>  ${logfile}
	#Pseudo-mapping with Kallisto to generate counts for RNAseq (from Pachter lab; eXpress/Cufflinks are disconstinued and not recommended/supported)
	command2="kallisto quant -t 20 -b 100 --single -l 200 -s 30 -i ${dir_index}/ITAG3.2_transcriptome -o ${dir_output} ${file} >> ${logfile} 2>&1"
	eval $command2 
	echo " *********** " >>  ${logfile}
	dir_output2=${dir_mapped2}/${base}
	mkdir -p ${dir_output2}
	ls -hl  ${file} >>  ${logfile}
	command3="STAR --readFilesCommand zcat --runThreadN 20 --genomeDir ${dir_genome} --readFilesIn ${file} --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignIntronMin 20 --alignIntronMax 10000 --bamRemoveDuplicatesType UniqueIdentical --outFilterMismatchNmax 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${dir_output2}/alignment >>  ${logfile} 2>&1"
	echo "$base STAR: mapping to tomato genome plus organelles"  >>  ${logfile}
	echo $command3  >>  ${logfile}
	eval $command3  
	echo " *********** " >>  ${logfile}
done