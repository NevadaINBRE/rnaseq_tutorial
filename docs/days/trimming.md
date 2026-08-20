

Following a QC analysis on sequencing results, one may detect stretches of low quality bases along reads, or a contamination by adapter sequence.
Depending on your research question and the software you use for mapping, you may have to remove these bad quality / spurious sequences out of your data.


**During this block, you will learn to :**

 * trim your data with fastp


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_03_trimming.pdf){target=_blank : .md-button }

[Fastp github](https://github.com/opengene/fastp){target=_blank : .md-button }


## to trim or not to trim ?

There are several ways to deal with poor quality bases or adapter contamination in reads, and several terms are used in the field, sometimes very loosely. We can talk about:

 * **Trimming**: to remove a part of, or the entirety of, a read (for quality reasons).
 	* **Hard trimming**: trim with a high threshold (eg. remove everything with QUAL<30).
 	* **Soft trimming**: trim with a low threshold (eg. remove everything with QUAL<10).
 * **Clipping**: to remove the end part of a read (typically because of adapter content).
 	* **Hard clipping**: actually removing the end of the read from the file (ie. with fastp).
 	* **Soft clipping**: ignoring the end of the read at mapping time (ie. what STAR does).


If the data will be used to perform **transcriptome assembly, or variant analysis, then it MUST be trimmed**.


In contrast, for applications based on **counting reads**, such as **Differential Expression analysis**, most aligners, such as [STAR](https://github.com/alexdobin/STAR), [HISAT2](http://daehwankimlab.github.io/hisat2/), [salmon](https://salmon.readthedocs.io/en/latest/salmon.html), and [kallisto](https://pachterlab.github.io/kallisto/manual), can handle bad quality sequences and adapter content by soft-clipping, and consequently they _usually_ do not need trimming.
In fact, **trimming can be detrimental** to the number of successfully quantified reads \[[William et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2)\].

Nevertheless, it is usually recommended to perform some amount of soft trimming (*eg.* [kallisto](https://www.biostars.org/p/389324/), [salmon](https://github.com/COMBINE-lab/salmon/issues/398) ).

If possible, we recommend to perform the mapping for both the raw data and the trimmed one, in order to compare the results for both, and choose the best.

**Question:** what could be a good metric to choose the best between the trimmed and untrimmed?

??? success "Answer"

	The number of uniquely mapped reads is generally what would matter in differential expression analysis. Of course, this means that you can only choose after you have mapped both the trimmed and the untrimmed reads.




## trimming with Fastp


The [fastp github](https://github.com/OpenGene/fastp/) gives very good examples of their software usage for both paired-end (`PE`) and single-end (`SE`) reads. We recommend you read their quick-start section attentively.



**Task 1:** 

 * Conduct a soft trimming on the mouseMT data

     - name the output folder : `030_d_trim/`.
     - unlike fastqc, you will have to launch fastp for each sample separately
     - fastp generates per-sample HTML and JSON QC reports automatically (no extra flags needed)


??? success "fastp sbatch script"

    We use fastp's default quality filtering, which performs a sliding-window quality cut from the right tail and auto-detects adapter sequences — no adapter file needed. Key options used:

     * **`-i`** : input fastq (single-end)
     * **`-o`** : trimmed output fastq
     * **`--thread 4`** : number of threads (match `--cpus-per-task`)
     * **`-j`** : JSON report (parsed by MultiQC)
     * **`-h`** : HTML report for manual inspection
     * **`--cut_right --cut_right_window_size 3 --cut_right_mean_quality 25`** : soft trim from 3' end when sliding-window mean quality drops below 25


    ```sh
    #!/usr/bin/bash
    #SBATCH --job-name=trim_mouseMT
    #SBATCH --time=01:00:00
    #SBATCH --cpus-per-task=4
    #SBATCH --mem=4G
    #SBATCH -o 030_l_trim_mouseMT.o
    #SBATCH --account=cpu-s5-biomarker_hunt-0
    #SBATCH --partition=cpu-core-0

    source ~/.bashrc
    conda activate rnaseq_env

    ## creating output folder, in case it does not exist
    mkdir -p 030_d_trim

    INPUT_FOLDER=/data/gpfs/assoc/biomarker_hunt/data/DATA/mouseMT

    fastp -i $INPUT_FOLDER/sample_a1.fastq \
          -o 030_d_trim/sample_a1.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_a1.json \
          -h 030_d_trim/030_l_trim_out.sample_a1.html

    fastp -i $INPUT_FOLDER/sample_a2.fastq \
          -o 030_d_trim/sample_a2.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_a2.json \
          -h 030_d_trim/030_l_trim_out.sample_a2.html

    fastp -i $INPUT_FOLDER/sample_a3.fastq \
          -o 030_d_trim/sample_a3.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_a3.json \
          -h 030_d_trim/030_l_trim_out.sample_a3.html

    fastp -i $INPUT_FOLDER/sample_a4.fastq \
          -o 030_d_trim/sample_a4.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_a4.json \
          -h 030_d_trim/030_l_trim_out.sample_a4.html

    fastp -i $INPUT_FOLDER/sample_b1.fastq \
          -o 030_d_trim/sample_b1.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_b1.json \
          -h 030_d_trim/030_l_trim_out.sample_b1.html

    fastp -i $INPUT_FOLDER/sample_b2.fastq \
          -o 030_d_trim/sample_b2.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_b2.json \
          -h 030_d_trim/030_l_trim_out.sample_b2.html

    fastp -i $INPUT_FOLDER/sample_b3.fastq \
          -o 030_d_trim/sample_b3.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_b3.json \
          -h 030_d_trim/030_l_trim_out.sample_b3.html

    fastp -i $INPUT_FOLDER/sample_b4.fastq \
          -o 030_d_trim/sample_b4.trimmed.fastq \
          --thread 4 \
          --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
          -j 030_d_trim/030_l_trim_out.sample_b4.json \
          -h 030_d_trim/030_l_trim_out.sample_b4.html
    ```

    On the cluster, this script is also in `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/030_s_trim.sh`


??? success "alternative sbatch using an array job"

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


	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=trim_mouseMT_array
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=4G
	#SBATCH -o 030_l_trim_mouseMT.%a.o
	#SBATCH --array=1-8%8
	#SBATCH --account=cpu-s5-biomarker_hunt-0
	#SBATCH --partition=cpu-core-0
	
	source ~/.bashrc
	conda activate rnaseq_env
	
	## creating output folder, in case it does not exist
	mkdir -p 030_d_trim
	
	INPUT_FOLDER=/data/gpfs/assoc/biomarker_hunt/data/DATA/mouseMT
	
	## each job grabs a specific line from sampleNames.txt
	SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sampleNames.txt)
	
	fastp -i $INPUT_FOLDER/${SAMPLE}.fastq \
	      -o 030_d_trim/${SAMPLE}.trimmed.fastq \
	      --thread 4 \
	      --cut_right --cut_right_window_size 3 --cut_right_mean_quality 25 \
	      -j 030_d_trim/030_l_trim_out.${SAMPLE}.json \
	      -h 030_d_trim/030_l_trim_out.${SAMPLE}.html
	```
	On the cluster, this script is also in	`/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/030bis_s_trim_array.sh`


**Task 2:** 

 * Use the the following script to run a QC analysis on your trimmmed reads and compare with the raw ones.


```sh
#!/usr/bin/bash
#SBATCH --job-name=trim-multiqc
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o 032_l_multiqc_trimmed.o
#SBATCH --account=cpu-s5-biomarker_hunt-0
#SBATCH --partition=cpu-core-0

source ~/.bashrc
conda activate rnaseq_env

## fastp already generated JSON/HTML QC reports; multiqc aggregates them directly
multiqc -n 032_r_multiqc_mouseMT_trimmed.html -f --title trimmed_fastp 030_d_trim/
```
On the cluster, you can find this script in : `/data/gpfs/assoc/biomarker_hunt/data/Solutions/mouseMT/032_s_multiqc_trimmed.sh`


