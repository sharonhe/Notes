# Download SRA data using command line utilities

   - Tools and instructions are available at: http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using

## Find the dataset to study
   - Go to http://www.ncbi.nlm.nih.gov/
   - Search "strawberry RNA-seq" in SRA database, select the dataset to study
   - Remember dataset number SRR1107997 in the following section: 
   ![image](https://cloud.githubusercontent.com/assets/16218822/11760246/c52364ce-a048-11e5-874a-3762197150c2.png)

## Download SRA files from NCBI using wget, and use fastq-dump to convert it to .fastq format
   - `wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR110/SRR1107997/SRR1107997.sra`
   - `nohup fastq-dump SRR1107997.sra`

## Download SRA files using fastq-dump in SRA Toolkit, the data will be automatically converted to .fastq format
   - `fastq-dump SRR390728`

# Usage of Samtools
## Installation of Samtools
- Go to http://www.htslib.org/download/ to download Samtools package
- Unzip it and open “INSTALL” file, according to the instruction, type:

```
cd samtools-1.2
make
make install # This avoids the $PATH setting step, and make the Samtools executable from any directory in terminal 

```
## Basic Samtools Commands
- Viewing all Samtools commands by simply typing:
`samtools`

- Calling a Samtool command, eg. flagstat, by typing:
`samtools flagstat`
It will show the Usage
For more information, type:
`samtools flagstat -help`

- Putting a “&” sign in the end of command allows the shell executes the command in the background. When finished, it will return a line:
`[2]-  Done   samtools sort Sample.bam Sample.sorted` 

### samtools view
adding “-h” allows the result display the header. The result will be automatically shown in .sam format, so this is also the usage of BAM-to-SAM conversion: 
`$ samtools view -h Sample.bam | more`
`$ samtools view -h Sample.bam > Sample.sam`

### SAM-to-BAM conversion
  - Checking the Usage of “SAM-to-BAM”:
`$ samtools view -help`
  - Execute “SAM-to-BAM” command: 
`samtools view -bT /Volumes/AdyBox/Shanshan_Bioinformatics/Reference_genome/hg_g1k_v37.fa Sample.sam > Sample.sam.bam`

### Samtools mpileup
- Only works for sorted and indexed bam files

- Usage: samtools mpileup [options] in1.bam [in2.bam [...]]

[Options] 
-f, —fasts-ref FILE    fax indexed reference sequence file
-v, --VCF            generate genotype likelihoods in VCF format
-u, —uncompressed    generate uncompressed VCF/BCF output
-g, —BCF             generate genotype likelihoods in BCF format

- Sample code:

```
$ samtools sort Exome_sample.bwamap.sam.bam Exome_sample.bwamap.sam.bam.sorted  # sort the bam file
$ samtools index Exome_sample.bwamap.sam.bam.sorted.bam  # index the bam file
$ samtools mpileup -f ../Reference_genome/hg_g1k_v37.fa Exome_sample.bwamap.sam.bam.sorted.bam > Exome_sample.mpileup   # generate genotype likelihoods in mpileup format
$ samtools mpileup -v -u -f ../Reference_genome/hg_g1k_v37.fa Exome_sample.bwamap.sam.bam.sorted.bam > Exome_sample.vcf   # generate genotype likelihoods in VCF format
$ samtools mpileup -g -f ../Reference_genome/hg_g1k_v37.fa Exome_sample.bwamap.sam.bam.sorted.bam > Exome_sample.bcf  # generate genotype likelihoods in BCF format

```
### samtools tview
- This is for viewing the reads aligned to the genome
- Usage: samtools view [options] <aln.bam> [ref.fasta]
[Options]:
   -d display      output as (H)tel or (C)uses or (T)ext 
   -p car:pos      go directly to this position
- Sample code:
`$ samtools tview -p 16:68194500 -d T Exome_sample.bwamap.sam.bam ../Reference_genome/hg_g1k_v37.fa | more`





# Download reference genome
## Download reference genome sequence (.fa) used in Galaxy
Instructions can be found at https://wiki.galaxyproject.org/Admin/UseGalaxyRsync. 
- List available genomes on server: 
`$ rsync --list-only rsync://datacache.g2.bx.psu.edu/indexes`
`$ rsync -- list-only rsync://datacache.g2.bx.psu.edu/indexes/hg19/`

- Download a file from the server to current directory. Note: using `sudo` and type password when necessary
`$ sudo rsync -avzP rsync://datacache.g2.bx.psu.edu/indexes/hg19/seq/hg19full.fa .`

## Download reference genome sequence (.gtf or .txt) from UCSC table browser
- Go to http://genome.ucsc.edu/cgi-bin/hgTables
- Put the information and get output
![image](https://cloud.githubusercontent.com/assets/16218822/11802141/0b117f92-a2a2-11e5-8a49-e66f81ae36cf.png)

# Usage of bedtools
## Installation of bedtools through homebrew
`$ brew install homebrew/science/bedtools`

## Commonly used bedtools function
- bedtools intersect
  - [Options] -wo: Write the original A and B entries plus the number of base pairs of overlap between the two features.
  - [Options] -split: each exon to be considered as separate interval
`$ bedtools intersect -wo -a GRCh38\:hg38.gtf -b Alus.bed.txt | more`
  - To check how many exons overlapped in these two files:
`$ bedtools intersect -wo -a GRCh38\:hg38.gtf -b Alus.bed.txt | wc -l`
or 
`$ bedtools intersect -split -wo -a GRCh38\:hg38.bed.txt -b Alus.bed.txt | wc -l`
  - To check how many genes overlapped in these two files:
`$ bedtools intersect -wo -a GRCh38\:hg38.gtf -b Alus.bed.txt | cut -f9 | cut -d “ “ -f2 | sort -u | wc -l` 
- bedtools bamtobed
  - [Options] -cigar: Add the CIGAR string to the BED entry as a 7th column
`$ bedtools bamtobed -cigar -i Sample.bam | more`
- bedtools bedtobam: 
  - Usage:   bedtools bedtobam [OPTIONS] -i <bed/gff/vcf> -g <genome>
  - Obtain the <genome> file used in this function through UCSC table browser:
    - Set the parameters in UCSC table browser as following:
![image](https://cloud.githubusercontent.com/assets/16218822/11825839/608c604c-a336-11e5-9ba9-a5e7939986cc.png)
    - Manipulate the genome file into the following format:
![image](https://cloud.githubusercontent.com/assets/16218822/11826035/ae91b386-a337-11e5-9a1b-a06a020f4c16.png)
  - Sample code: 
`$ bedtools bedtobam -i hg38.gtf -g ../Shanshan_Bioinformatics/Reference_genome/hg38.genome > hg38.bam`
- bedtools getfasta
  - Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>
  - Sample code:
    - view result per gene
`$ bedtools getfasta -fi ../Shanshan_Bioinformatics/Reference_genome/hg19full.fa -bed RefSeq_bed.txt -fo Refseq_Chr22_hg19.gtf.fasta`
    - view result per feature (eg. exon)
`$ bedtools getfasta -split -fi ../Shanshan_Bioinformatics/Reference_genome/hg19full.fa -bed RefSeq_bed.txt -fo Refseq_Chr22_hg19.gtf.fasta`

# Usage of bowtie2
## Installation of bowtie2
1. Go to http://bowtie-bio.sourceforge.net/bowtie2/faq.shtml to download the latest release
2. Add Bowtie 2 directory to my PATH environment variable (See the markdown file “ToolsSetup.md”)
## bowtie2-build
- bowtie2-build builds a Bowtie index from a set of DNA sequences.
- Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
- Sample code:
`$ bowtie2-build hg19full.fa hg19`
## bowtie2 for mapping
### Usage: 
bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

  <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gzip’ed (extension: .gz) or bzip2’ed (extension: .bz2).
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip’ed (extension: .gz) or bzip2’ed (extension: .bz2).
  <r>        Files with unpaired reads.
             Could be gzip’ed (extension: .gz) or bzip2’ed (extension: .bz2).
  <sam>      File for SAM output (default: stdout)
Performance:
  -p/—threads <int> number of alignment threads to launch

### Sample code:
-end-to-end mapping:
`$ bowtie2 -p 4 -x ../Reference_genome/hg19_build/hg19 ChIP_Sample.fastq -S ChIP_Sample.sam`

-local mapping:
`$ bowtie2 -p 4 --local -x ../Reference_genome/hg19_build/hg19 ChIP_Sample.fastq -S ChIP_Sample_local.sam`

# Usage of BWA
## Installation of BWA
1. Go to http://sourceforge.net/projects/bio-bwa/files/ to download the latest release
2. Untar the file, cd to bwa folder in terminal, and type `make`
3. Add bwa directory to my PATH environment variable (See the markdown file “ToolsSetup.md”)

## bwa index
To use BWA, you need to first index the reference genome
`$ bwa index hg_g1k_v37.fa`

## bwa mem
### Usage: 
bwa mem [options] <idxbase> <in1.fq> [in2.fq]
[Options] -t INT: number of threads
### Sample code for mapping pair-end reads:
`$ bwa mem -t 4 ../Reference_genome/hg_g1k_v37.fa Exome_sample_R1.fastq Exome_sample_R2.fastq > Exome_sample.bwamap.sam`

# Usage of bcftools
## Installation of bcftools
1. Go to http://www.htslib.org/download/ to download the latest release
2. Untar the file, cd to bcftools folder in terminal, and type `make install`
## bcftools call
### About bcftools call   
SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.

### Usage:
bcftools call [options] <in.vcf.gz>
[Options]
-v, —variants-only             output variant sites only
-m, —multiallelic-caller       alternative model for multi allelic and rare-variant calling (conflicts with -c)
-O, —output-type <b|u|z|v>     output type: ‘b’ compressed BCF; ‘u’ uncompressed BCF; ‘z’ compressed VCF; ‘v’ uncompressed VCF [v]
-o, —output <file>             write output to a file [standard output]

### Sample code

```

$ bcftools call -v -m -O z -o Exome_sample.vcf.gz Exome_sample.vcf  # Variant calling
$ gunzip -c Exome_sample.vcf.gz # for viewing the above output file, zcat won’t work on OS X.
$ gunzip -c Exome_sample.vcf.gz | grep -v “^#” | wc -l   # To count how many variants are there

```


