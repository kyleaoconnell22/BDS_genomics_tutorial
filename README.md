# Perform whole genome alignment and variant calling

---------------------------------
## Overview of Page Contents

+ [Purpose](#P)
+ [Dependencies](#D)
+ [Workflow Overview](#WFO)
+ [Workflow Steps](#WFS)
+ [Input Files](#IN)
+ [Download Files](#DOW)
+ [Tutorial](#TUT)
+ [Running on GCP](#GCP)
+ [License](#LIC)
+ [Contacts](#CON)

## **Purpose** <a name="P"></a>

This workflow is designed to benchmark genomic applications on cloud platforms. This example is run on GCP with CPUs, but for the tutorial we will run everything locally. This tools requires Python 3. You can check your python version by typing (in the terminal) python --version
If you have Python 2, you can create a conda env (see install info below) for python 3.

## **Dependencies** <a name="D"></a>

This workflow requires a number of dependencies to be installed on your local machine or cloud VM. For the purpose of this tutorial, the easiest way to install everything will be with the Anaconda distribution.

### Install miniconda package

### *Installing conda on a Mac:* 
The curl command is used to download the installer from the web. Note that the -O flag is a capital o not a zero.

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

Install miniconda into $HOME/miniconda3  
* Type 'yes' to agree to the license
* Press Enter to use the default install directory
* Type 'yes' to initialize the conda install

`bash Miniconda3-latest-MacOSX-x86_64.sh`

Refresh your terminal session to see conda

`bash`

Test that conda is installed. Will print info about your conda install.

`conda info`


### *Installing conda on a Linux:*  
Fetch the miniconda installer with wget  
`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

Install miniconda into $HOME/miniconda3  
* Type 'yes' to agree to the license
* Press Enter to use the default install directory
* Type 'yes' to initialize the conda install  

`bash Miniconda3-latest-Linux-x86_64.sh`

Refresh your terminal session to see conda.  
`bash`

Test that conda is installed. Will print info about your conda install.  
`conda info`

### *Creating python 3 environment if needed:*  
`conda create --name py3 python=3.5`  
`conda init`  
`source activate py3`    
   
## **Workflow Overview** <a name="WFO"></a>

The workflow is run by the script execute_mapping.py which can be run locally, or via sh startup script on GCP.

### *Pipeline Dependencies*
+ **bwa**:  [bwa](http://bio-bwa.sourceforge.net/) is used for genome alignment. 
+ **GATK**: [GATK](https://gatk.broadinstitute.org/hc/en-us) a collection of tools from the Broad Institute used throughout the workflow.
+ **Samtools**: [Samtools](http://www.htslib.org/) is used for manipulating and analyzing SAM and BAM files, could be used parallel to most gatk tools. Here we use it to calculate read depth and other metrics
+ **VCFtools**: [VCFtools](http://vcftools.sourceforge.net/) is used for manipulating and analyzing VCF files
+ **docker**: used for running [deepvariant](https://github.com/google/deepvariant). It is a pretty big package so if you don't want to install it locally it isn't essential for this tutorial. It will already be on GCP (try 'conda install -c conda-forge docker', otherwise just install from 'https://docs.docker.com/engine/install/ubuntu/')
+ **numpy**: python package needed to run a few of the scripts.  

   *Conda install commands*  
   `conda install -c bioconda bwa`  
   `conda install -c bioconda gatk4`  
   `conda install -c bioconda samtools=1.11`  
   `conda install -c bioconda vcftools`  
   `conda install -c conda-forge numpy` 


### *Mandatory Arguments*  
+ **path_to_data**: Absolute path to FASTQ data (input>subject>files). FASTQ data should be paired end sequence files (R1 and R2). See tutorial for example of dir structure
+ **path_to_ref**: Absolute path for reference genome. This does not have to be indexed, but it saves a lot of time if it is. In this tutorial the reference has already been indexed.
+ **path_to_results**: Absolute path for results (default is current directory)
+ **path_to_resources**: Absolute path for GATK resources directory, used for BQSR and VQSR
+ **threads**: Number of threads
+ **variant_caller**: Variant caller can be either haplotypecaller or deepvariant. Haplotype caller is from the GATK package. Deepvariant is Google's deep learning variant caller
+ **filter_vcf**: In the case of haplotypecaller, what tool would you like to use to filter variants? Options include VQSR, CNN (in progress), or NONE

### *Optional Arguments*  
+ **index_genome**: Option to use bwa index to index the reference genome. This step is time consuming so we will download an indexed genome
+ **verbose**: Print debugging information
+ **deleteResultFiles**: Delete results directory after execution

## **Input Files** <a name="IN"></a>

+ **Fastq data**: These are the files outputted from the sequencer. They are usually filtered to remove Illumina adapters and contamination. These need to be both forward and reverse reads for this pipeline to work. The files we will use in the tutorial have already been filtered.  
+ **Reference genome**: Here we are using GRCh38 human reference, which is the most recent iteration of the human reference genome. Most of this workflow can be run with or without alt chromosomes, but for VQSR to work you have to use without alt chromsomes because all the resource VCFs were generated with the GATK Resources Reference and the site names need to match up. Alt chromosomes are additional versions of the chromosomes included in the main reference (which can only by definition have two for each chromosome since humans are diploid). For example, some human populations may all have a large insertion that is not included in the main reference, and this would be captured in the alt chromosome sequence, but would be lacking from the main reference at that site.
+ **Resources files**: We will download the GATK resource files used for BQSR and VQSR filtering. 


## **Workflow Steps** <a name="WFS"></a>

![Workflow](https://github.com/kyleaoconnell22/BDS_genomics_tutorial/blob/main/workflow_figures/Genome_alignment_workflow.png)


### **BWA mem**

Align filtered fastq files to a reference genome. The script is designed to move through several subjects, and each subject can have multiple fastq file pairs. For example, if you have the directory 'input', it could contain subject1, subject2, and subject3, and each of these could have fastq1.R1.fastq (forward reads), fastq1.R2.fastq (reverse reads) and then fastq2.R1.fastq, fastq2.R2.fastq etc. (The reason there may be multiple fastq files for each subject (in this case person) is that deep sequencing (30-50x or more) may require more sequencing than any one machine can output. Another issues is that PCR-free libraries require several libraries to be prepared for each subject to have enough DNA to sequence deeply. Spreading the sequence across a few libraries and a few sequencing machines also distributes error (mainly levels of duplication) more widely, and should give better sequence quality overall).

This step also adds read group labels, which define groups of reads generated by a single run of a sequencing instrument. The details are beyond the scope of this tutorial, but basically bwa is going to add a tag with barcode, sequencing instrument, and sample name to the header of each read in the SAM/BAM file. This info is used later when we merge BAMs from different sequencing runs for the same subject to mark duplicates and recalibrate quality scores for each sequencing run. In order to get the read group label right, the naming of the fastq files is fairly specific for the script as currently written. See the tutorial. 

Also note that BWA can take a long time for deeply sequenced subjects, but in the toy data used in the tutorial it runs super fast. For example, the 50x data we are using has 71M reads, and takes about 20 hours on GCP with 32 threads. 

### **GATK SortSam**

BWA will output a SAM file, which is a file with information about reads mapped to a reference (BAM is binary equivalent, CRAM are column-oriented binary container equivalent). Sorting is where you sort by the order of the reads in the reference genome (for example, Chromosome 1, then Chromosome 2 etc.). This is necessary for future steps and speeds up analyses.

If you have multiple fastq pairs for a subject, then bwa mem and SortSam will be run on each pair separately and then we merge the fastq pairs for the same subject in the next step.

### **GATK MarkDuplicates** 

This step is doing two things. First, it is bringing together all the different SAM files and merging them into a single BAM file (but maintaining the read group info for each read so you know which fastq file it came from originally). Second, it marks reads as duplicates, using the read group tags to identify which reads come from which sequencing run. A duplicate tag is added to the header so that these reads can be excluded for variant calling.  

Read duplicates can arise via two processes. First, duplicates can come from PCR duplication when the libraries are generated in the lab (library or PCR duplicates). The data we are using here is PCR-free, so that is not a concern in this case. PCR-free means that the data was generated using several input libraries then sequenced several times. The second kind of duplication arises due to sequencing error, where a single amplification cluster is detected incorrectly as multiple clusters. These are called optical duplicates. Optical duplicates can by definition only exist from a single sequencing run which is why we need the read group tags to differentiate fastq files from multiple runs. 

Also note that optical duplicates are only an issue with Illumina 2500 machines, newer sequencing machines (now the NovaSeq which is where our data comes from here) use a different kind of flowcell that has lower incidence of duplication. 

### **Decide on variant caller**
At this point in the script, things change depending on if you are using GATK haplotypecaller or Google deepvariant. If haplotypecaller, then you will conduct base quality score recalibration, then variant calling and variant filtering. If deepvariant, then you first index your bam file then call variants directly after MarkDuplicates with no additional filtering.

### **Base Quality Score Recalibration**
For a really good explanation of what is happening here, go read the [Broad Institute description](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) . In a nutshell, sequencing machines assign quality scores to how confident the machine is about the identity of a given base. These estimates have systematic error, and GATK applies a machine learning model to try to correct for some of this error. 

There are two commands used to do this, gatk BaseRecalibrator and gatk ApplyBQSR. BaseRecalibrator builds a model for the 'corrected' scores, and ApplyBQSR modifies the BAM file according to these adjusted scores. Note that this part of the pipeline also uses the read group tags, and models quality scores for each read group separately.

### **Call germline variants**

Germline variants are those that are inherited through the germline (sperm and egg). Somatic variants are those that arise in an individual (e.g. cancer-associated variants) that are not passed on to the next generation. Somatic variant calling is a whole different pipeline that we will not discuss further here. 

### **GATK haplotypecaller**
Haplotypecaller is a germline snp (single nucleotide polymorphism) and indel (insertion/deletion) caller. It works by comparing the BAM file to the reference genome. Anywhere it identifies variation, it reassembles (think realign) the sequence for that region and then calls variants where the subject BAM differs from the reference sequence. This difference could be one bp (SNP) or several bp missing from the reference (insertion) or missing from the subject (deletion).  
At this point, we expect about 3M variants for a deeply-sequenced human genome. There may be many more though due to false positives. This step has high sensitivity (most true variants included), but low specificity (lots of false positives included too).  

### ***GATK Variant Call Score Recalibration***
Similar to BQSR, here we are applying a machine learning model to filter variants, but instead of using an existing model built from known patterns from general Illumina runs (BQSR), this step builds a new model based on paramaters modeled from example data that you feed the program. For example, if you have known sites from existing human genome projects, VQSR builds a model that describes these sites, then rescores your variants based on how well they match the parameter space of the known data the model was built from.  
Similar to BQSR above, this analysis has two steps, VariantRecalibrator, which builds the model, and 
Apply VQSR, which rescores your VCF file according to the model. Note that it changes the Pass/Fail flag from Pass to LowVQSLOD (low score), so the number of variants in the VCF file does not change, and these have to be hard filtered at the very end of any pipeline. 
For a more detailed description, go to the [GATK man page](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-)

### ***DeepVariant***
DeepVariant uses a CNN-based program that produces pileup image tensors from each BAM file, then classifies each tensor using the CNN model, and finally calls variants from the pileups. 

This Google [blog post](https://google.github.io/deepvariant/posts/2020-02-20-looking-through-deepvariants-eyes/) does a great job of explaining the method in more detail.  Although no additional machine learning filtering is needed as with haplotypecaller, there will be a lot of low quality variants in the deepvariant VCF that need to be hard filtered.

## **Download Files** <a name="DOW"></a>  

### ***Download Repository and CD to base dir***
`git clone https://github.com/Deloitte/genomics_tutorial.git`  
`cd genomics_tutorial-main`  
For those in a hurry, all directories and files can be setup and downloaded with the helper script. Otherwise, you can go the long route and do everything step by step below. 

`sh scripts/setup_tutorial.sh`

Congrats! Assuming everything went well, skip ahead to the Tutorial section below. If you had any issues then consider just running the rest of the Downloads section manually.

### ***Download Reference Genome***

But first, set up directory structure if you are taking the long route  

`mkdir resources reference input results input/HG00096 results/HG00096 output_examples`

And then move the input data to the correct place

`mv input.fastq.zip input/HG0096 | unzip input/HG00096/input.fastq.zip`


All resources we download will come from the [GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false). We are going to just download one resource today for the tutorial, but in the full pipeline we end up using all of these. 

`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict -output reference/Homo_sapiens_assembly38.dict`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta -output reference/Homo_sapiens_assembly38.fasta`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt -output reference/Homo_sapiens_assembly38.fasta.64.alt`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb -output reference/Homo_sapiens_assembly38.fasta.64.amb`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann -output reference/Homo_sapiens_assembly38.fasta.64.ann`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt -output reference/Homo_sapiens_assembly38.fasta.64.bwt`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac -output reference/Homo_sapiens_assembly38.fasta.64.pac`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa -output reference/Homo_sapiens_assembly38.fasta.64.sa`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai -output reference/Homo_sapiens_assembly38.fasta.fai`  

From GCP bucket:  
`gsutil -m cp  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"  
  reference`  

### *Download Resources*

From URL:  

`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz`  
`curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi`  

From GCP bucket:  
`gsutil -m cp -r   
  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"  
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"   


## **Tutorial** <a name="TUT"></a>

Make sure you download all the files and set up the directory structure first! If you have not done that then go back up one level and get everything set up. For those in a hurry, you can also just run:  
'sh setup.sh'

### *Run Commands*  

You could run the whole pipeline by running the python script like this:  
`python3 scripts/execute_mapping.py input reference/Homo_sapiens_assembly38.fasta --path_to_resources resources --variant_caller haplotypecaller --filter_vcf VQSR --threads 4`

But, in this tutorial we are going to run the commands individually so that you can get a feel for what is happening at each step:

### Assign global variables
`RESULTS='results/HG00096'`  
`RESOURCES='resources'`  
`REF='reference/Homo_sapiens_assembly38.fasta'`  
`EXOUT='output_examples'`  

### Map fastq files to reference genome with BWA mem to produce a SAM file with mapping info
`bwa mem -t 2 -R '@RG\tID:HFHJKDSXX.L004\tLB:lib1\tPL:LLUMINA\tPU:HFHJKDSXX.L004.CGGACAAC_TCCGGATT\tSM:HG00096' reference/Homo_sapiens_assembly38.fasta input/HG00096/HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L004_001.R1.fastq.gz input/HG00096/HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L004_001.R2.fastq.gz > $RESULTS/out.sam`

### View the Samfile, note the chromosome order  
`samtools view $RESULTS/out.sam | head -20`

### SortSam and convert to BAM  
`gatk SortSam -I $RESULTS/out.sam -O $RESULTS/HG00096.bam -SO coordinate`

### View sorted BAM, note chromosome order now  
`samtools view $RESULTS/HG00096.bam | head -20`

### Mark Duplicates  
`gatk MarkDuplicates I=$RESULTS/HG00096.bam O=$results $RESULTS/HG00096.markdup.bam M=$RESULTS/metrics.txt`

### View output, just look at top rows to see % duplication rate split by PCR and optical duplicates, here it is very low as expected in PCR-free libraries  
`head $RESULTS/metrics.txt`  

### Index the BAM, this is only needed for deepvariant  
`gatk BuildBamIndex -I $RESULTS/HG00096.markdup.bam -O $RESULTS/HG00096.markdup.bai`  

### *GATK haplotypecaller Steps*

### First calculate model for recalibrating base quality scores  
`gatk BaseRecalibrator -I $RESULTS/HG00096.markdup.bam -O $RESULTS/bqsr_out.txt --known-sites $RESOURCES/Homo_sapiens_assembly38.known_indels.vcf.gz -R $REF`

### Look at output of model, notice our read group info is used to sort by read group. These are the new scores for each base  
`cat $RESULTS/bqsr_out.txt`

### Apply BQSR model  
`gatk ApplyBQSR -R $REF -I $RESULTS/HG00096.markdup.bam -O $RESULTS/HG00096.bam  -bqsr-recal-file $RESULTS/bqsr_out.txt`

### Lets see how well mapping went and our read depth going into variant calling  
#### Samtools flagstat will give us the % of reads that mapped  
`samtools flagstat $RESULTS/HG00096.bam`  

#### Samtools depth + python script will give us the read depth  
`samtools depth $RESULTS/HG00096.bam > $RESULTS/depth.out`
`python scripts/find_mean.py $RESULTS/depth.out`  

### Call snps and indels with haplotype caller, this will take about 6min on a local computer, so if you are pressed for time an example output is in the output_examples dir
`gatk HaplotypeCaller -I $RESULTS/HG00096.bam -O $RESULTS/HG00096_haplotypecaller.vcf -R $REF --native-pair-hmm-threads 4 --output-mode EMIT_VARIANTS_ONLY`

### Filter raw SNPs with VQSR, to save time we are going to ignore indels here, but in a real analysis you would run it for snps and then indels sequentially.  
### First build model, here we are using a VCF of the whole genome from the output_examples, because with the toy dataset there are not enough sites to actually build the model
### Run locally on my laptop took about ~30 min. 
`gatk VariantRecalibrator -V $EXOUT/HG00096_haplotypecaller.vcf --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an FS -an MQ -an SOR -an DP -mode SNP -resource:hapmap,known=false,training=true,truth=true,prior=15 resources/hapmap_3.3.hg38.vcf.gz -resource:omni,known=false,training=true,truth=true,prior=12 resources/1000G_omni2.5.hg38.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10 resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=7 resources/Homo_sapiens_assembly38.dbsnp138.vcf -O $RESULTS/snps.recal --tranches-file $RESULTS/snps.tranches`

### View output files  
`head -100 output_examples/HG00096/snps.recal`  
`cat output_examples/HG00096/snps.tranches`  

### Apply model, this only takes 1.3 min on my local machine
`gatk ApplyVQSR -R $REF -V $EXOUT/HG00096_haplotypecaller.vcf -O $RESULTS/HG00096_filt.vcf --tranches-file $RESULTS/snps.tranches --recal-file $RESULTS/snps.recal -mode SNP --exclude-filtered  -ts-filter-level 90.0`  

### See how many variants are in your toy data vcf file from haplotype caller  
`vcftools --vcf $EXOUT/HG00096_haplotypecaller.vcf`

*Notice we have one individual, and in our subsampled data the model found 123 SNPs/Indels*

In the full dataset we would have ~5M variants, but after filtering we have 3.6-4.2M.

### Call variants using deepvariant
This runs three steps, the longest of which is training (23 min on my local machine), then call variants (15min) and then postprocess variants (5min)

`BIN_VERSION="1.1.0"`  
`S='HG00096'`  
`INPUT_DIR="${PWD}/results/${S}"`  
`REF_DIR="${PWD}/reference"`  
`REF="Homo_sapiens_assembly38.fasta"`  
`BAM="${S}.markdup.bam"`  
`OUTPUT_DIR="${INPUT_DIR}"`  
`OUTPUT_VCF="${S}_deepvariant.vcf"`  

`sudo docker run -v "${INPUT_DIR}":"/input" -v "${OUTPUT_DIR}":"/output" -v "${REF_DIR}":"/refdir" google/deepvariant:"1.1.0" /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref="/refdir/${REF}"  --reads="/input/${BAM}"  --output_vcf="/output/${OUTPUT_VCF}"  
 --num_shards=4  
    --intermediate_results_dir /output/intermediate_results_dir`  

### Check our results  
`vcftools --vcf $EXOUT/${S}_deepvariant.vcf`

Open and explore the [deepvariant visual report](https://github.com/kyleaoconnell22/BDS_genomics_tutorial/blob/main/output_examples/HG00096_deepvariant.visual_report_example.html)

## **Running on GCP** <a name="GCP"></a>
1) Create a balanced persistent disk to hold the data and command line programs. For example, you could install miniconda here and use that to install other programs, or you can install them one by one. 
This disk will need at least 600GB of storage for a single 30x subject. Scale up according to sequencing depth and number of subjects. We used debian-10-buster-v20210316 for machine type, but this can be flexible.
This VM will need to have all dependencies installed. You will also need a directory with input fastqs, and a directory with the reference genome, as well as a resources dir with files (except the reference which needs to be in the reference directory) from the [GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)

2) Create gs buckets for the following categories
+ **scripts**: This bucket needs to have the python script "execute_mapping.py", and also the DeepVariant sh script "run-deepvariant.sh". Both are included in the [scripts](https://github.com/kyleaoconnell22/BDS_genomics_tutorial/tree/main/scripts) directory
+ **startup**: This bucket needs to have "launch_mapping.sh" An example is included in the [scripts](https://github.com/kyleaoconnell22/BDS_genomics_tutorial/tree/main/scripts) directory
+ **output_logs***: Once the pipeline finishes it will write the log file to here based on the experiment name in your gcloud command
+ **outputVCFs**: Once the pipeline finishes it will write the output VCFs here for each subject, and in the case of DeepVariant, will output a nice html summary of variants called.

3) Launch the script with a command like this:
+ **Sub out $string for your specific strings**
+ **We are using ec2-standard-32 for machine type, so it has 32 cores and 128 GB RAM. You could scale this down for testing, but for 30x genomes and up you need the threads and the RAM**

gcloud compute --project=$ProjectName instances create $InstanceName  --zone=us-east4-c --machine-type=e2-standard-32 --subnet=usgcpadvnpd2 --no-address --maintenance-policy=MIGRATE --service-account=933976732552-compute@developer.gserviceaccount.com --scopes=https://www.googleapis.com/auth/devstorage.read_write,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append --image=$ImageName --image-project=$ProjectName --boot-disk-size=$DiskSize --boot-disk-type=pd-balanced --boot-disk-device-name=instance-2 --no-shielded-secure-boot --shielded-vtpm --shielded-integrity-monitoring --reservation-affinity=any --metadata experimentname=ExampleExpName2021,mappingthreads=32,variantcaller=haplotypecaller,filter=VQSR,startup-script-url=gs://bds-genomics-test-bucket/startup/launch_mapping.sh

## **License** <a name="LIC"></a>

GNU Lesser General Public License v3.0


## **Contacts:** <a name="CON"></a>

Kyle O'Connell, PhD
kyoconnell@deloitte.com

Collin Lobb, PhD
clobb@deloitte.com

