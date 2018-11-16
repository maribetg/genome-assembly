# genome-assembly
steps for genome assembly


# emerge two files for velvet
./RAD/velvet_1.2.10/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl ./Genome-Siphonoperla-torrentium/FCH5H32CCXY_L4_wHADPI063284-113_1.fq ./Genome-Siphonoperla-torrentium/FCH5H32CCXY_L4_wHADPI063284-113_2.fq genome.fq

./RAD/velvet_1.2.10/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl ./Genome-Siphonoperla-torrentium/genome-second-run/FCHLJHLCCXY_L8_WHRDINSqgmRAAAAA-701501_1.fq ./Genome-Siphonoperla-torrentium/genome-second-run/FCHLJHLCCXY_L8_WHRDINSqgmRAAAAA-701501_2.fq genome02.fq

# Run velveth
velveth velvelhgenome 31,45,2 -fastq -shortPaired genome.fq 
velvetg velvelhgenome_31 -exp_cov auto 

velveth velvelhgenome 31,40,2 -fastq -shortPaired genome02.fq 
velvetg velvelhgenome_31 -exp_cov auto 

velveth velvelhgenome 31,40,2 -fastq -shortPaired ./Genome-Siphonoperla-torrentium/genome-clean-alldata.fastq 
velvetg velvelhgenome_31 -exp_cov auto 

velveth velvelhgenome 81 -fastq -shortPaired ./Genome-Siphonoperla-torrentium/genome-clean-alldata.fastq 
velvetg velvelhgenome_81 -exp_cov auto 

# installing sga
./sga/src/configure --with-sparsehash=/sparsehash --with-bamtools=/usr/local --prefix=/sga/sga/install
make

# sga for quality control
./sga/SGA/sga preprocess --pe-mode 1 FCH5H32CCXY_L4_wHADPI063284-113_1.fq FCH5H32CCXY_L4_wHADPI063284-113_2.fq > genome.fastq

./sga/SGA/sga preprocess --pe-mode 1 -q -f FCH5H32CCXY_L4_wHADPI063284-113_1.fq FCH5H32CCXY_L4_wHADPI063284-113_2.fq > genome-clean.fastq

./sga/SGA/sga preprocess --pe-mode 1 --quality-trim=40 --quality-filter=3 --permute-ambiguous --remove-adapter-fwd=ACACTCTTTCCCTACACGACGCTCTTCCGATCT --remove-adapter-rev=AGACGTGTGCTCTTCCGATCT FCH5H32CCXY_L4_wHADPI063284-113_1.fq FCH5H32CCXY_L4_wHADPI063284-113_2.fq > ./sga/genome-clean02.fastq

./sga/SGA/sga preprocess --pe-mode 1 --quality-trim=40 --quality-filter=3 --permute-ambiguous --remove-adapter-fwd=ACACTCTTTCCCTACACGACGCTCTTCCGATCT --remove-adapter-rev=AGACGTGTGCTCTTCCGATCT ./genome-second-run/FCHLJHLCCXY_L8_WHRDINSqgmRAAAAA-701501_1.fq ./genome-second-run/FCHLJHLCCXY_L8_WHRDINSqgmRAAAAA-701501_2.fq > ./sga/genome-clean-second-run.fastq

./sga/SGA/sga preprocess --pe-mode 1 --quality-trim=40 --quality-filter=3 --permute-ambiguous --remove-adapter-fwd=ACACTCTTTCCCTACACGACGCTCTTCCGATCT --remove-adapter-rev=AGACGTGTGCTCTTCCGATCT ./genome-second-run/FCHLJHLCCXY_L8_WHRDINSqgmRAAAAA-701501_1.fq ./genome-second-run/FCHLJHLCCXY_L8_WHRDINSqgmRAAAAA-701501_2.fq ./genome-first-run/FCH5H32CCXY_L4_wHADPI063284-113_1.fq ./genome-first-run/FCH5H32CCXY_L4_wHADPI063284-113_2.fq> genome-clean-alldata.fastq

./SGA/sga index -a ropebwt -t 3 genome-clean02.fastq
./SGA/sga index -a ropebwt -t 3 genome.fastq
./SGA/sga index -a ropebwt -t 3 genome-clean-alldata.fastq

./SGA/sga preqc -t 8 genome-clean02.fastq > genome-clean02.preqc
./SGA/sga preqc -t 8 genome.fastq > genome.preqc
./sga/SGA/sga preqc -t 8 genome-clean-alldata.fastq > genome-clean-alldata.preqc

python2 ./sga/src/bin/sga-preqc-report.py genome-clean02.preqc ./sga/src/examples/preqc/*.preqc 
python2 ./sga/src/bin/sga-preqc-report.py genome.preqc ./sga/src/examples/preqc/*.preqc
python2 ./sga/src/bin/sga-preqc-report.py genome-clean-alldata.preqc ./sga/src/examples/preqc/*.preqc

# kmer size
./kmergenie/kmergenie-1.7048/kmergenie ./sga/genome-clean02.fastq -o ./kmergenie/genome-clean02
./kmergenie/kmergenie-1.7048/kmergenie genome-clean-alldata.fastq -o ./kmergenie/genome-clean-alldata #K was run between 15 to 121
./kmergenie/kmergenie-1.7048/kmergenie genome-clean-alldata.fastq -l 10 -o ./kmergenie/genome-clean-alldata10 # change minimum k to 10

# assembly with minia
./minia-v2.0.7-bin-Linux/bin/minia -in ./sga/genome-clean02.fastq -kmer-size 19 -abundance-min 22 -out ./minia/minia-genome-clean02-19
./minia-v2.0.7-bin-Linux/bin/minia -in ./sga/genome-clean02.fastq -kmer-size 21 -abundance-min 21 -out ./minia/minia-genome-clean02-21
./minia-v2.0.7-bin-Linux/bin/minia -in genome-clean-alldata.fastq -kmer-size 15 -abundance-min 21 -out ./minia/minia-genome-clean-alldata-15



## Alignment mitogenomes and genome
mkdir mapping
#convert fastq to fasta. Using fastz-toolkit
fastq_to_fasta -i genome.fq -o genome.fasta

# normalize fasta file
picard-tools NormalizeFasta I=./minia/assembly-41/minia-assem01.contigs.fa O=./minia/assembly-41/minia-assem01.contigs.nor.fa

1.Creating a reference directory using bowtie2
bowtie2-build ./minia/assembly-41/minia-assem01.contigs.nor.fa ./mapping/bowtie-index
bowtie2 -x bowtie-index -q stonefly-mitgen-fastq.fq -S align-mito-local.sam --local
bowtie2 -x bowtie-index -f stonefly-mitgen-seq-genebank.fasta -N 1 -S align-bowtie-fasta.sam

#t ransform sam output of bowtie to .bam
./samtools-1.8/samtools view -b -h -S ./mapping/align-mito-local.sam -o ./mapping/align-mito-local.bam

# clean up some reads
./samtools-1.8/samtools fixmate -O bam ./mapping/align-mito-local.sam ./mapping/align-mito-local-fix.bam

# sort .bam files
./samtools-1.8/samtools sort ./mapping/align-mito-local-fix.bam -o ./mapping/align-mito-local-fix-sorted.bam

# index .bam files
./samtools-1.8/samtools index ./mapping/align-mito-local-fix-sorted.bam

# view alignment
./samtools-1.8/samtools tview ./mapping/align-mito-local-fix-sorted.bam


## Create a index file for bwa
bwa index ./minia/assembly-41/minia-assem01.contigs.nor.fa -p ./mapping/bwa-index
bwa mem bwa-index stonefly-mitgen-fastq.fq > align-bwa.sam
./samtools-1.8/samtools view -b -h -S ./mapping/align-bwa.sam -o ./mapping/align-bwa.bam
./samtools-1.8/samtools sort ./mapping/align-bwa.bam -o ./mapping/align-bwa-sorted.bam
./samtools-1.8/samtools index ./mapping/align-bwa-sorted.bam
./samtools-1.8/samtools tview ./mapping/align-bwa-sorted.bam

## Seeing alignments in igv
type igv in consola

## Pot using circos
