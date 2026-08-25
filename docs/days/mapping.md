
Once you are happy with your read sequences in your FASTQ files, you can use a mapper software to align the reads to the genome and thereby find where they originated from.


**At the end of this lesson, you will be able to :**

 * identify the differences between a local aligner and a pseudo aligner.
 * perform genome indexing appropriate to your data.
 * map your RNA-seq data onto a genome.



## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_04_Mapping.pdf){target=_blank : .md-button }

[STAR website](https://github.com/alexdobin/STAR){target=_blank : .md-button }



## Building a reference genome index

Before any mapping can be achieved, you must first *index* the genome want to map to. 

To do this with STAR, you need two files:

 * a *fasta* file containing the sequences of the chromosome (or genome contigs)
 * a *gtf* file containing annotations (ie. where the genes and exons are)

We will be using the Ensembl references, with their accompanying GTF annotations.

!!! note

	While the data are already on the server here, in practice or if you are following this course without a teacher,
	you can grab the reference genome data from the [Ensembl ftp website](https://www.ensembl.org/info/data/ftp/index.html).

	In particular, you will want a mouse [DNA fasta file](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/) and [gtf file](http://ftp.ensembl.org/pub/current_gtf/mus_musculus/).

  Take note of the genome sequence and annotation versions, you will need this in your paper's methods section!


**Task :** Using STAR, build a genome index for the mouse mitochondrial chromosome.

 * .fasta and .gtf files are in : `/data/gpfs/assoc/biomarker_hunt/data/DATA/Mouse_MT_genome/`.
   * we will need to create a MT specific GTF file 
 * create the index in the folder `041_d_STAR_mouseMT_reference`
 * the aligner is available in the `rnaseq_env` conda environment.
 * this job should require less than 4Gb and 10min to run. 

!!! info "STAR basic parameter for genome index generation"

	From the [manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf). Refer to it for more details

	 * `--runMode genomeGenerate` : running STAR in index generation mode
	 * `--genomeDir </path/to/genomeDir>` : output folder for the index
	 * `--genomeFastaFiles </path/to/genome/fasta1>` : chromosome sequences fasta file (can be several files)
	 * `--sjdbGTFfile </path/to/annotations.gtf>` : annotation gtf file
	 * `--runThreadN <NumberOfThreads>` : number of threads to run on 
	 * `--sjdbOverhang <ReadLength-1>` : length of the genomic sequence around the annotated junctions to be used in constructing the splice junctions database. Ideally : read length - 1.

	 Additionally, because the genome is so small here (we only use the mitochondrial chromosome after all), you will need the following advanced option:

	 * `--genomeSAindexNbases 5` : must be scaled to `min(14, log2(GenomeLength)/2 - 1)`, so 5 in our case


!!! note
	
	While your indexing job is running, you can read ahead in STAR's manual to prepare the next step : mapping your reads onto the indexed reference genome.


!!! example "practical"

    Use the command line utility `awk` to create a copy of the GTF file in your `mouseMT` folder filtered for **MT** GTF entries. 
    Since `awk` is a lightweight program, we can run this filtering step on the headnode

??? success "GTF filter"

	```sh
	G_GTF=/data/gpfs/assoc/biomarker_hunt/data/DATA/Mouse_MT_genome/Mus_musculus.GRCm39.116.gtf

        awk '/^#/ || $1 == "MT"' $G_GTF > mt_only.gtf
	```


??? success "STAR indexing script"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-build
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=3G
	#SBATCH -o 041_l_star_index.o
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	G_FASTA=/data/gpfs/assoc/biomarker_hunt/data/DATA/Mouse_MT_genome/Mus_musculus.GRCm39.dna.chromosome.MT.fa
	G_GTF=mt_only.gtf

	mkdir -p 041_d_STAR_mouseMT_reference

	STAR --runMode genomeGenerate \
	     --genomeDir 041_d_STAR_mouseMT_reference \
	     --genomeFastaFiles $G_FASTA \
	     --sjdbGTFfile $G_GTF \
	     --runThreadN 4 \
	     --genomeSAindexNbases 5 \
	     --sjdbOverhang 99

	```

	It can be found on the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/041_s_star_index.sh`

**Extra task :** Determine how you would add an additional feature to your reference, for example for a novel transcript not described by the standard reference.

??? success "Answer"

	You would edit the gtf file to add your additional feature(s), following the [proper format](https://www.ensembl.org/info/website/upload/gff.html).



!!! note "Note"

	In case you've got multiple FASTA files for your genome (eg, 1 per chromosome), you may just list them with the `genomeFastaFiles` option as follow:

	`--genomeFastaFiles /path/to/genome/fasta1.fa /path/to/genome/fasta2.fa /path/to/genome/fasta3.fa ...`


## Mapping reads onto the reference


**Task :** Using STAR, align the raw FASTQ files of the mouseMT dataset against the mouse mitochondrial reference you just created.

 * if were not able to complete the previous task, you can use the index in `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/041_d_STAR_mouseMT_reference` .
 * search the STAR manual for the option to output a BAM file sorted by coordinate.
 * search the STAR manual for the option to output a geneCounts file.
 * put the results in folder `042_d_STAR_map_raw/` .


!!! info "STAR basic parameters for mapping"

	Taken again from the manual:

	- `--genomeDir </path/to/genomeDir>` : folder where you have put the genome index
	- `--readFilesIn </path/to/read1> ` : path to a fastq file. If the reads are paired, then also include the path to the second fastq file
	- `--runThreadN <NumberOfThreads>`: number of threads to run on.
	- `--outFileNamePrefix  <prefix> ` : prefix of the output files, typically something like `output_directory/sampleName` . 



!!! Note

	Take the time to read the parts of the [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) which concern you: a bit of planning ahead can save you a lot of time-consuming/headache-inducing trial-and-error on your script.


!!! Warning

	Mapping reads and generating a sorted BAM from one of the mouseMT FASTQ file will take less than a minute and very little RAM, but on a real dataset it should take from 15 minutes to an hour per sample and require at least 30GB of RAM.



??? success "STAR mapping script"

	We will be using a job array to map each file in different job that will run at the same time.

	First create a file named `sampleNames.txt`, containing the sample names:

	```
	sample_a1
	sample_a2
	sample_a3
	sample_a4
	sample_b1
	sample_b2
	sample_b3
	sample_b4
	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/sampleNames.txt`

	Then for our script:

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=00:10:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=1G
	#SBATCH -o 042_l_STAR_map_raw.%a.o
	#SBATCH --array 1-8%8
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	mkdir -p 042_d_STAR_map_raw

	SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sampleNames.txt)

	FASTQ_NAME=/data/gpfs/assoc/biomarker_hunt/data/DATA/mouseMT/${SAMPLE}.fastq

	STAR --runThreadN 4 --genomeDir 041_d_STAR_mouseMT_reference \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix 042_d_STAR_map_raw/${SAMPLE}. \
	     --quantMode GeneCounts \
	     --readFilesIn $FASTQ_NAME

	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/042_s_STAR_map_raw.sh`

	and its results can be found at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/042_d_STAR_map_raw/`


	The options of STAR are :

	 * `--runThreadN 4 ` : 4 threads to go faster.
	 * `--genomeDir 041_STAR_reference` : path of the genome to map to.
     * `--outSAMtype BAM SortedByCoordinate ` : output a coordinate-sorted BAM file.
     * `--outFileNamePrefix 042_STAR_map_raw/${SAMPLE}.` : prefix of output files.
     * `--quantMode GeneCounts` : will create a file with counts of reads per gene.
     * `--readFilesIn $FASTQ_NAME` : input read file.



## QC on the aligned reads

You can call MultiQC on the STAR output folder to gather a report on the individual alignments.


**Task :** use `multiqc` to generate a QC report on the results of your mapping.

 * Evaluate the alignment statistics. Do you consider this to be a good alignment?
 * How many unmapped reads are there? Where might this come from, and how would you determine this?
 * What could you say about library strandedness ? 

??? success "script and answers"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=map-multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o 043_l_multiqc_map_raw.o
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	multiqc -n 043_r_multiqc_mouseMT_mapped_raw.html -f --title mapped_raw 042_d_STAR_map_raw/
	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/043_s_multiqc_map_raw.sh`


	[ Download the report ](../assets/html/043_multiqc_mouseMT_mapped_raw.html){target=_blank : .md-button }


## Comparison of mapping the trimmed reads

After having mapped the raw reads, we also map the trimmed reads and then compare the results to decide which one we want to use for the rest of our analysis.


[ trimmed reads mapping  report ](../assets/html/045_multiqc_mouseMT_mapped_trimmed.html){target=_blank : .md-button }



??? success "scripts for the mapping of trimmed reads"
	
	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=00:10:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=1G
	#SBATCH -o 044_l_STAR_map_trimmed.%a.o
	#SBATCH --array 1-8%8
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	mkdir -p 044_d_STAR_map_trimmed

	SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sampleNames.txt)

	FASTQ_NAME=030_d_trim/${SAMPLE}.trimmed.fastq

	STAR --runThreadN 4 --genomeDir 041_d_STAR_mouseMT_reference \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix 044_d_STAR_map_trimmed/${SAMPLE}_trimmed. \
	     --quantMode GeneCounts \
	     --readFilesIn $FASTQ_NAME
	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/044_s_STAR_map_trimmed.sh`


	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=map-trim-multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o 045_l_multiqc_mouseMT_mapped_trimmed.o
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	multiqc -n 045_r_multiqc_mouseMT_mapped_trimmed.html -f --title mapped_trimmed 044_d_STAR_map_trimmed/
	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/045_s_multiqc_mouseMT_mapped_trimmed.sh`


## QC report of mapping for the Liu2015 and Ruhland2016 dataset

**Liu2015**

Take the time to look at the following reports: 

[ Liu2015 raw reads mapping  report ](../assets/html/042_r_star_aln_raw_QC_Liu2015.html){target=_blank : .md-button }
[ Liu2015 trimmed reads mapping  report ](../assets/html/044_r_star_map_trim_QC_Liu2015.html){target=_blank : .md-button }

Which one would you choose? 


**Ruhland**

[ Ruhland2016 raw reads mapping  report ](../assets/html/034_r_STAR_multiqc_Ruhland2016.html){target=_blank : .md-button }

??? success "scripts for Ruhland2016"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-build
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=10G
	#SBATCH -o 041_l_star_index.o
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	G_FASTA=/data/gpfs/assoc/biomarker_hunt/data/DATA/Mouse_MT_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa
	G_GTF=/data/gpfs/assoc/biomarker_hunt/data/DATA/Mouse_MT_genome/Mus_musculus.GRCm39.116.gtf

	mkdir -p 041_d_STAR_Ruhland_reference

	STAR --runMode genomeGenerate \
	     --genomeDir 041_d_STAR_Ruhland_reference \
	     --genomeFastaFiles $G_FASTA \
	     --sjdbGTFfile $G_GTF \
	     --runThreadN 4 \

	```
	
	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=00:10:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=8G
	#SBATCH -o 044_l_STAR_map_trimmed.%a.o
	#SBATCH --array 1-8%8
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	mkdir -p 044_d_STAR_map_trimmed

	SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sampleNames.txt)

	FASTQ_NAME=030_d_trim/${SAMPLE}.trimmed.fastq.gz

	STAR --runThreadN 4 --genomeDir 041_d_STAR_Ruhland_reference \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix 044_d_STAR_map_trimmed/${SAMPLE}_trimmed. \
	     --quantMode GeneCounts \
	     --readFilesIn $FASTQ_NAME
	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/Ruhland2016/044_s_STAR_map_trimmed.sh`


	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=map-trim-multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o 045_l_multiqc_ruhland2016_mapped_trimmed.o
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	multiqc -n 045_r_multiqc_ruhland_mapped_trimmed.html -f --title mapped_trimmed 044_d_STAR_map_trimmed/
	```
	it can also be found in the cluster at `/data/gpfs/assoc/biomarker_hunt/data/Solutions/Ruhland2016/045_s_multiqc_mouseMT_mapped_trimmed.sh`
