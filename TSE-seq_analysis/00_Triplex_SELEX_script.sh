#!/bin/sh

#  Run as: ./Triplex_SELEX_Script.sh
#  Created by Sara Wernig-Zorc on 05.02.2021
#  5 samples: 2x input TFO library and 3x Triplex RNA

#OUTPUT COLORS
GREEN="\033[1;32m"
PURPLE="\033[1;35m"
BLUE="\033[1;34m"
RESET="\033[1;0m"

set -e

echo "$PURPLE--------- Triplex-SELEX analysis ---------$RESET"

echo "$PURPLE--------- Set project related variables ---------$RESET"

project_folder=/Volumes/Sara_RAID_20TB/Sara_data_2018/Triplex-SELEX_project
samples="$project_folder/raw_data/*_1.fastq.gz"

echo  "$PURPLE The project folder: $RESET" $project_folder
echo  "$PURPLE Annotation file: $RESET" $annotation
echo  "$PURPLE Genome: $RESET" $genome
echo  "$PURPLE Genome fasta file: $RESET" $fasta
echo  "$PURPLE Samples: $RESET" $(basename ${samples/%_1.fastq.gz/})

cd $project_folder

echo  "$PURPLE--------- START OF ANALYSIS ---------$RESET"

mkdir -p raw_data processed_data/trimmed processed_data/aligned logs scripts results/fastQC

echo  "$GREEN--------- Quality Control ---------$RESET"

fastqc -t 10 -o $project_folder/results/fastQC $project_folder/raw_data/*.fastq.gz
multiqc $project_folder/results/fastQC/ -o $project_folder/results/fastQC/

echo  "$BLUE--------- Adapter trimming and read pairing ---------$RESET"

cd $project_folder/raw_data/
####2.1 Cutadapt
for file in *R1_001.fastq.gz; do echo "Processing file: " $file; cutadapt --trim-n -a GATCGGAAGAGCACACGTCTG -a AGAGCACACGTCTG $file | cutadapt -u 3 -a A{10}N{90} --no-indels -e 0.16666666666666666 - | cutadapt -O 8 --match-read-wildcards -g GTTCAGAGTTCTACAGTCCGACGATCSSS -m 18 -o ${file/%.fastq.gz/_trimmed.fastq.gz} - ; done

for file in *R2_001.fastq.gz; do echo "Processing file: " $file; cutadapt --trim-n --match-read-wildcards -n 2 -g T{100} -a SSSGATCGTCGG -m 18 -o ${file/%.fastq.gz/_trimmed.fastq.gz} $file - ; done

for file in *fastq.gz; do echo "File: " $file; gzip -dc $file ; done

#2.2.1 Trimmomatic
for file in *R1_001_trimmed.fastq; do java -jar /Users/sara/tools/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 12 -phred33 -summary ${file/%.fastq/_step2.log} $file ${file/%.fastq/_cropped.fastq} ILLUMINACLIP:/Volumes/Sara_RAID_20TB/Sara_data_2018/Triplex-SELEX_project/logs/NGS_contaminants.fa:5:10:10 TRAILING:10 MINLEN:20 CROP:27; done

for file in *R2_001_trimmed.fastq; do java -jar /Users/sara/tools/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 12 -phred33 -summary ${file/%.fastq/_step2.log} $file ${file/%.fastq/_cropped.fastq} ILLUMINACLIP:$project_folder/logs/NGS_contaminants.fa:5:10:10 MINLEN:20 CROP:27; done

#2.2.2 PE Trimmomatic
#for forward in *R1_001_trimmed.fastq ; do java -jar /Users/sara/tools/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -threads 20 $forward ${forward/%R1_001_trimmed.fastq/R2_001_trimmed.fastq} ${forward/%R1_001_trimmed.fastq/_01_paired.fastq} ${forward/%R1_001_trimmed.fastq/_01_unpaired.fastq} ${forward/%R1_001_trimmed.fastq/_02_paired.fastq} ${forward/%R1_001_trimmed.fastq/_02_unpaired.fastq} ILLUMINACLIP:/Volumes/Sara_RAID_20TB/Sara_data_2018/Triplex-SELEX_project/logs/NGS_contaminants.fa:5:10:10  MINLEN:20 CROP:27; done

mv *trimmed* $project_folder/processed_data/trimmed/
mv *log $project_folder/logs/

echo  "$GREEN--------- Extract lines that begin with CCTCCTTTTTTCTTTTTTTT ---------$RESET"

for file in $project_folder/processed_data/trimmed/*trimmed_cropped.fastq; do echo $file; grep -B1 -A2 "^CCTCCTTTTTTCTTTTTTTT" $file | grep -v "^--$" - > ${file/%.fastq/_CCTCCT.fastq}; done


echo  "$GREEN--------- Number of reads == 27nt ---------$RESET"

cd $project_folder/processed_data/trimmed/
for file in *_CCTCCT.fastq; do echo $file; java -jar /Users/sara/tools/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 12 -phred33 -summary ${file/%.fastq/_step3.log} $file ${file/%.fastq/_27nt.fastq} MINLEN:27 ; done


echo  "$GREEN--------- Extract only the variable region ---------$RESET"

#For MEME motif enrichment
for file in *_CCTCCT_27nt.fastq; do echo $file; ~/tools/fastX_toolkit/bin/fastx_trimmer -i $file -o ${file/%_27nt.fastq/_varReg.fastq} -f 21 ; done


echo  "$GREEN--------- Fastq to Fasta format ---------$RESET"

#for file in *.gz; do echo $file;  gzip -d $file; done

for file in $project_folder/processed_data/trimmed/*_CCTCCT_27nt.fastq; do echo $file; gsed -n '1~4s/^@/>/p;2~4p' $file > ${file/%.fastq/.fasta}; done

echo  "$GREEN--------- Merge replicates ---------$RESET"

cat input*.fasta > input-TFO-library-merged_trimmed_cropped_CCTCCT_varReg.fasta
cat triplex*.fasta > triplex-RNA-merged_trimmed_cropped_CCTCCT_varReg.fasta

echo  "$GREEN--------- Read count after cleaning ---------$RESET"

for file in *_CCTCCT_27nt.fastq; do echo "File: " $file;  awk 'NR%4==2{c++; l+=length($0)} END{print "Number of reads: "c; print "Number of bases in reads: "l }' ; done

echo  "$GREEN--------- MEME: Motif enrichment ---------$RESET"

/Users/sara/tools/MEME_v5.3.2/bin/meme $project_folder/processed_data/varRegion/triplex-RNA-rep01_S3_L001_R1_001_trimmed_cropped_CCTCCT_varReg.fasta -neg $project_folder/processed_data/varRegion/input-TFO-library-rep01_S1_L001_R1_001_trimmed_cropped_CCTCCT_varReg.fasta -objfun de -minw 7 -rna -nmotifs 10 -o $project_folder/results/MEME/rep01/

/Users/sara/tools/MEME_v5.3.2/bin/meme $project_folder/processed_data/varRegion/triplex-RNA-rep02_S4_L001_R1_001_trimmed_cropped_CCTCCT_varReg.fasta -neg $project_folder/processed_data/varRegion/input-TFO-library-rep01_S1_L001_R1_001_trimmed_cropped_CCTCCT_varReg.fasta -objfun de -minw 7 -rna -nmotifs 10 -o $project_folder/results/MEME/rep02/

/Users/sara/tools/MEME_v5.3.2/bin/meme $project_folder/processed_data/varRegion/triplex-RNA-rep03_S5_L001_R1_001_trimmed_cropped_CCTCCT_varReg.fasta -neg $project_folder/processed_data/varRegion/input-TFO-library-rep02_S2_L001_R1_001_trimmed_cropped_CCTCCT_varReg.fasta -objfun de -minw 7 -rna -nmotifs 10 -o $project_folder/results/MEME/rep03/

#Merged replicates

/Users/sara/tools/MEME_v5.3.2/bin/meme $project_folder/processed_data/varRegion/triplex-RNA-merged_trimmed_cropped_CCTCCT_varReg.fasta -neg $project_folder/processed_data/varRegion/input-TFO-library-merged_trimmed_cropped_CCTCCT_varReg.fasta -objfun de -minw 7 -rna -nmotifs 10 -o $project_folder/results/MEME/merged

echo  "$GREEN--------- MAIN ANALYSIS ---------$RESET" 

#1.
#for file in *R1_001.fastq.gz; do echo "Processing file: " $file; cutadapt -u 3  $file | cutadapt -q 10 -o ${file/%.fastq.gz/_trimmed.fastq.gz} - ; done

#(optional)
#for file in *fastq.gz; do echo "File: " $file; gzip -dc $file ; done

#2.
#for file in *_trimmed.fastq; do echo $file; gsed -n '1~4s/^@/>/p;2~4p' $file > ${file/%.fastq/.fasta}; done

#3.
$project_folder/scripts/01_Extract_triplex.R

#4.
$project_folder/scripts/02_Triplex-SELEX.Rmd

echo  "$GREEN--------- SALMON: Read quantification - PART I ---------$RESET"

#Create in silico transcriptome
salmon index -i RNA_TFO_library -k 7 -p 20 -t $project_folder/processed_data/in_silico_library/Random_DNA_seq_Nr.1e+05.fasta

#Quantify reads Control
salmon quant -i RNA_TFO_library/ -l SF -r $project_folder/raw_data/FASTA/input-TFO-library-rep01+rep02_trimmed_filt20.fasta -o $project_folder/processed_data/readQuant/ctrl/ --noFragLengthDist --noSingleFragProb  --noEffectiveLengthCorrection --noLengthCorrection --maxOccsPerHit 100000 --allowDovetail --threads 10

#Quantify reads Triplex
salmon quant -i RNA_TFO_library/ -l SF -r $project_folder/raw_data/FASTA/triplex-RNA-rep02+rep03_trimmed_filt20.fasta -o $project_folder/processed_data/readQuant/triplex/ --noFragLengthDist --noSingleFragProb  --noEffectiveLengthCorrection --noLengthCorrection --maxOccsPerHit 100000 --allowDovetail --threads 10

echo  "$GREEN--------- SALMON: Read quantification - PART II ---------$RESET"

#In silico Library
awk 'NR%2{printf "%s ",$0;next;}1' $project_folder/processed_data/in_silico_library/Random_DNA_seq_Nr.1e+05.fasta > $project_folder/processed_data/trimmed/RNA_inSilico_library-oneLine.csv

gsed s'/ /  /'g $project_folder/processed_data/trimmed/RNA_inSilico_library-oneLine.csv | gsed s'/>//'g > $project_folder/processed_data/trimmed/RNA_inSilico_library-Sequence.csv

mkdir $project_folder/results/SALMON/

#Add sequences to input lib
join -1 1 -2 1 <(sort -k1,1 $project_folder/processed_data/readQuant/ctrl/quant.sf) <(sort -k1,1 $project_folder/processed_data/trimmed/RNA_inSilico_library-Sequence.csv) > $project_folder/results/SALMON/input_quant_wSeq.csv

#Add sequences to triplex lib
join -1 1 -2 1 <(sort -k1,1 $project_folder/processed_data/readQuant/triplex/quant.sf) <(sort -k1,1 $project_folder/processed_data/trimmed/RNA_inSilico_library-Sequence.csv) > $project_folder/results/SALMON/triplex_quant_wSeq.csv

#logFC

gsed 's/\t/\n/g' top_enriched_triplexes.txt | gsed 's/Random/>Random/g' - > top_enriched_triplexes.fasta

echo  "$BLUE--------- DONE ---------$RESET"
