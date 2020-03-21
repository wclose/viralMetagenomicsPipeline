## Metagenomic Community Analysis of Virus-like Particles

This repo can be used to generate all of the desired output files for community analysis (shared OTUs, taxonomies, alpha/beta diversity, ordiations, etc.). The workflow is designed to work with [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Conda](https://docs.conda.io/en/latest/) with minimal intervention by the user. The workflow is designed to be run on a high-performance computing cluster (HPC) because some of the steps are very computationally intense and require signficant resources (>100 GB of RAM and >24 cores). If you do not have access to a cluster, you should be able to adjust the settings to make it execute locally though it will admittedly be difficult. An overview of the workflow is provided in the [rulegraph](rulegraph.svg).

This workflow was designed to analyze viral community composition using DNA isolated from virus-like particles in mouse feces. It has been used to successfully analyze data from 380 samples with 10 Gb of sequencing per sample on an Illumina NovaSeq but should be able to easily scale up or down depending on the amount of data being used. Samples from other environments should also work fine though you will need to update the contaminate indices used for [metagenomeBowtie2Dependencies.sh](code/bash/metagenomeBowtie2Dependencies.sh) and [metagenomeDecontaminateReads.sh](code/bash/metagenomeDecontaminateReads.sh) execution.

### Overview

The steps below are meant to serve as a brief overview of the pipeline. For more information, please see the comments throughout the various scripts and the documentation for the tools used.

**1. Read processing**
* **Purpose**: Remove low-quality sequences and any potential host contaminants (mouse and human in this case).
* Trim read pairs and remove low-quality regions ([Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)).
* Remove potential contaminant sequences ([Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)).

<br />

**2. Assembly**
* **Purpose**: Combine quality-controlled reads into contigs.
* Independently assemble reads from each sample into contigs ([metaSPAdes](http://cab.spbu.ru/software/spades/)).
* Combine all generated contigs into a single contig library for the entire dataset ([Bash](https://www.gnu.org/software/bash/manual/bash.html)).
* Remove duplicate contigs and containments to reduce redundancy of sequences ([BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)).
* Remove contigs below a specific size threshold to improve downstream processing ([SeqKit](https://github.com/shenwei356/seqkit)).

<br />

**3. Viral contig curation**
* **Purpose**: Remove low-quality contigs and any contigs predicted as being non-viral.
* Predict which contigs are viral ([VirSorter](https://github.com/simroux/VirSorter) and [VirFinder](https://github.com/jessieren/VirFinder)).
* Remove any non-viral contigs from the contig library ([Bash](https://www.gnu.org/software/bash/manual/bash.html)).

<br />

**4. Binning**
* **Purpose**: Group contigs together into discrete units based on similarity to each other and distance from others.
* Map reads from each sample to the viral contig library to generate coverage and abundance statistics ([BWA-MEM](http://bio-bwa.sourceforge.net/)).
* Combine contigs into metagenomic bins based on tetranucleotide frequencies, coverage, and abundance correlations ([MetaBat2](https://bitbucket.org/berkeleylab/metabat/src/master/)).
* Remove any bins that contain non-viral marker genes ([CheckM](https://github.com/Ecogenomics/CheckM/wiki)).
* Combine all bins together into a single viral bin library ([Bash](https://www.gnu.org/software/bash/manual/bash.html)).

<br />

**5. Shared OTU file creation**
* **Purpose**: Summarize the number of counts for each viral bin across each sample in the dataset. Bins are treated as operational taxonomic units (OTUs).
* Map reads from each sample to the viral bin library ([BWA-MEM](http://bio-bwa.sourceforge.net/)).
* Calculate number of reads per bin for each sample. ([Samtools](http://www.htslib.org/)).
* Combine all read counts into a single shared OTU file (R - [tidyverse](https://www.tidyverse.org/)).
* Rarefy the shared OTU files by subsampling to a specified threshold (R - [tidyverse](https://www.tidyverse.org/)).
* Normalize shared OTU counts based on the total length of each bin (R - [tidyverse](https://www.tidyverse.org/)).

<br />

**6. Taxonomic classification**
* **Purpose**: Provide basic taxonomic information for each of the viral OTUs (bins).
* Classify metagenome assembled bins (MAGs; [CAT/BAT](https://github.com/dutilh/CAT)).

<br />

**7. Calculate diversity metrics**
* **Purpose**: Calculate various metrics needed for comparing community compositions/properties.
* Calculate alpha (within sample) and beta (between sample) diversity metrics (R - [tidyverse](https://www.tidyverse.org/), [furrr](https://github.com/DavisVaughan/furrr), [vegan](https://cran.r-project.org/web/packages/vegan/index.html)).
* Generate principal coordinates of analysis (PCoA) and non-metric multidimensional (NMDS) ordinations using beta diversity (R - [tidyverse](https://www.tidyverse.org/), [vegan](https://cran.r-project.org/web/packages/vegan/index.html), [ecodist](https://cran.r-project.org/web/packages/ecodist/ecodist.pdf)).

<br />

**8. Quality control**
* **Purpose**: Identify any stages where read or contig qualities are inadequate.
* Map reads to contaminant sequences to determine degree of contamination ([FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)).
* Summarize read quality scores, adapter content, etc. ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
* Summarize contig assembly quality scores, lengths, etc. ([Quast](http://quast.sourceforge.net/)).
* Combine quality control diagnostics together into a single interactive report ([MultiQC](https://multiqc.info/)).

### Usage

#### Dependencies
* MacOSX or Linux operating system.
* Have access to a high-performance computing (HPC) cluster. 
* Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and configure for use on a cluster. Instructions for installing `conda` on a cluster are provided by me [HERE](https://github.com/um-dang/conda_on_the_cluster).
* Have paired end sequencing data (fastq.gz).
> **NOTE:** The workflow assumes you are using the Swift Accel-NGS 1S Plus DNA Library Kit (Swift Biosciences; cat. no. 10024) for library prepartion and sequencing on an Illumina machine. If using a different library kit, you may need to adjust the `metagenomeHeadcrop` setting in [config.yaml](config/config.yaml). If using a different sequencer, you will need to adjust the adapter types/location during execution of [metagenomeTrimmomatic.sh](code/bash/metagenomeTrimmomatic.sh) in the [Snakefile](Snakefile).

<br />

#### Running analysis

**1.** Clone this repository and move into the project directory.
```
git clone https://github.com/wclose/viralMetagenomicsPipeline.git
cd viralMetagenomicsPipeline
```

<br />

**2.** Transfer all of your raw paired-end sequencing data into `data/raw/`. 
> **NOTE:** Because of the way sample names are automatically parsed, files are assumed to be in the form `140142_TATGCCAG-TTCCATTG_S154_R1_001.fastq.gz`. Everything before the first underscore is assumed to be the sample name (Ex: `140142`). Files are also assumed to be gzipped fastq files ending in `fastq.gz`. You will need to adjust accordingly if this is not the case.

<br />

Copy sequences to the raw data directory.
```
cp PATH/TO/SEQUENCEDIR/* data/raw/
```

<br />

**3.** Create the master `snakemake` environment.
> **NOTE:** If you already have a conda environment with snakemake installed, you can skip this step.
```
conda env create -f envs/yaml/snakemake.yaml
```

<br />

**4.** Activate the environment that contains `snakemake`.
```
conda activate snakemake
```

<br />

**5.** Edit the options in [config.yaml](config/config.yaml) file to set downstream analysis options.
```
nano config/config.yaml
```

<br />

Things to change (everything else can/should be left as is):
* **metagenomeControls**: Put the names (just the names of the samples, not the full filename with all of the sequencer information) of all of your controls here. Make sure names don't have underscores as these will not be parsed correctly.
* **metagenomeHeadcrop**: The number of bases to trim from each read with Trimmomatic. This setting depends on the library preparation kit used and your sequencing conditions. You may need to increase/decrease this depending on the FastQC results in the MultiQC reports generated as part of the pipeline. 
* **metagenomeVirfinderFpr** and **metagenomeVirfinderFdr**: The false-positive and false-discovery rates to be used for filtering VirFinder results when predicting viral contigs. 
* **metagenomeSubthresh**: The number of reads to subsample to for each sample when generating the shared OTU count file.
* **metagenomeAlpha**: The desired alpha (within sample) diversity metrics to calculate.
* **metagenomeBeta**: The desired beta (between sample) diversity metrics to calculate.

<br />

**6.** Test the workflow to make sure everything looks good. This will print a list of the commands to be executed without running them and will error if something isn't set properly.
```
snakemake -np
```

<br />

**7.** If you want to see how everything fits together, you can run the following to generate a flowchart of the various steps. Because of the amount of parallel processes, I highly recommend using the rulegraph as opposed to the directed acyclic graph (DAG) output from `snakemake` as the DAG will be extremely cluttered. I have included the pre-generated [rulegraph](rulegraph.svg) for the pipeline within this repository. You may need to download the resulting image locally to view it properly. 
> **NOTE:** If you are using MacOSX, you will need to install `dot` by installing `graphviz` with [homebrew](https://brew.sh/) or some alternative process before running the following command.
```
snakemake --rulegraph | dot -Tsvg > rulegraph.svg
```

<br />

**8.** Before running any jobs on the cluster, change the `ACCOUNT` and `EMAIL` fields in the following files. If you need to change resource requests (increase/decrease time, memory, etc.), you can find those settings in the cluster profile configuration files as well.
* Slurm: [Cluster profile configuration](config/slurm/cluster.yaml) and the [cluster submission script](code/slurm/snakemake.sh).

<br /> 

**9.** Run the Snakemake workflow. Snakemake will submit each task as an individual job on the cluster, monitor the jobs, and submit new jobs as needed without needing any intervention from the user. Should there be an error, Snakemake will tell you where the error occurred and halt all jobs. It is recommended to submit the workflow using the [cluster submission script](code/slurm/snakemake.sh) as this will create a job that manages the workflow for you. 
```
sbatch code/snakemake.sh
```

<br />

**10.** Assuming that no errors occur, the workflow will create multiple directories within the `data/` directory. Outputs of interest will be as follows:
* **`data/shared/`**: This directory will contain output files containing normalized counts of each metagenomic bin (viral OTU) across all samples.
* **`data/diversity/`**: This directory contains all metrics for measuring diversity including the alpha diversity table and any beta diversity distance files.
* **`data/catbat/`**: This directory contains the taxonomic assignments of the various metagenomic bins (viral OTUs) using [CAT/BAT](https://github.com/dutilh/CAT) and the [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) viral database.
* **`data/multiqc/`**: This directory contains all of the qualtiy control reports with interactive graphs. Reports are generated at different checkpoints throughout the pipeline and I highly encourage you to look at them, adjust the pipeline settings as necessary, and rerun the workflow to generate the best quality data.
