
Read counting refers to the quantification of an *"expression level"*, or abundance, from reads mapped onto a reference genome/transcriptome.
This *expression level* can take several forms, such as a count, or a fraction (RPKM/FPKM/TPM), and concern different entities (exons, transcripts, genes) depending on your biological application.


**During this lesson, you will learn to:**

 * differentiate between different levels of counting and their relevance for different questions.


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_05_ReadCounting.pdf){: .md-button }

[featureCounts website](http://subread.sourceforge.net/featureCounts.html){: .md-button }


## Read counting with featureCounts

<!-- Suggestion: Perhaps add a quick note on why use featureCounts when we've already quantified using STAR? -->

The [featureCount website](http://subread.sourceforge.net/featureCounts.html) provides several useful command-line examples to get started.
For more details on the algorithm behavior (with multi/overlapping reads for instance), you can refer to the package's [User's guide](http://subread.sourceforge.net/SubreadUsersGuide.pdf) (go to the read summarization chapter).


**Task :** 

 * Decide which parameters are appropriate for counting reads from the mouseMT dataset. Assume you are interested in determining which genes are differentially expressed.
 * Count the reads from your trimmed BAM files using featureCounts.
 * How do the featureCount-derived counts compare to those from STAR ?
 * featureCount requirements : 400M RAM / BAM file
 * featureCount requirements : 2 min CPU time / BAM file

??? done "featureCounts raw script"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=featurecount
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=4G
	#SBATCH -o 050_l_featureCounts_mouseMT.o
	#SBATCH -e 050_l_featureCounts_mouseMT.e
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	G_GTF=mt_only.gtf

	inFOLDER=042_d_STAR_map_raw

	mkdir -p 050_d_featureCounts_raw_mouseMT

	featureCounts -T 8 -a $G_GTF -t exon -g gene_id \
		-o 050_d_featureCounts_raw_mouseMT/050_r_featureCounts_raw_mouseMT.counts.txt \
		$inFOLDER/*.bam

	```


??? done "featureCounts trimmed script"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=featurecount
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=4G
	#SBATCH -o 050_l_featureCounts_mouseMT.o
	#SBATCH -e 050_l_featureCounts_mouseMT.e
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	G_GTF=mt_only.gtf

	inFOLDER=044_d_STAR_map_trimmed

	mkdir -p 051_d_featureCounts_trimmed_mouseMT

	featureCounts -T 8 -a $G_GTF -t exon -g gene_id \
		-o 051_d_featureCounts_trimmed_mouseMT/050_r_featureCounts_trimmed_mouseMT.counts.txt \
		$inFOLDER/*.bam

	```



Review the output matrix. Is this informative? What are the informative columns? Could we get more information? These identifiers are great, but what about gene information?

Let's view the information in the GTF file and see how we can add this information in this output matrix.

Now rerun featureCounts with the `--extraAttributes` option configured with the new output filename as `052_d_featureCounts_trimmed_extraAttributes_mouseMT/051_r_featureCounts_mouseMT.counts.extraAttributes.txt`.


??? done "featureCounts trimmed with extraAttributes script"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=featurecount
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=4G
	#SBATCH -o 051_l_featureCounts_extraAttributes_mouseMT.o
	#SBATCH -e 051_l_featureCounts_extraAttributes_mouseMT.e
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0

	source ~/.bashrc
	conda activate rnaseq_env

	G_GTF=mt_only.gtf

	inFOLDER=044_d_STAR_map_trimmed

	mkdir -p 052_d_featureCounts_trimmed_extraAttributes_mouseMT

	featureCounts -T 8 -a $G_GTF -t exon -g gene_id \
	    --extraAttributes "transcript_id,gene_name,gene_source,gene_biotype,transcript_name,transcript_biotype,tag" \
		-o 052_d_featureCounts_trimmed_extraAttributes_mouseMT/051_r_featureCounts_trimmed_mouseMT.counts.extraAttributes.txt \
		$inFOLDER/*.bam

	```
