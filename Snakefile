# Snakefile
# William L. Close
# Schloss Lab
# University of Michigan

# Purpose: Snakemake workflow for analyzing viral metagenomic sequences purified from virus-like particles

# Path to config file
configfile: "config/config.yaml"

# Function for aggregating list of raw sequencing files.
metagenomeSamples = list(set(glob_wildcards(os.path.join('data/raw/', '{sample}_{readNum, R[12]}_001.fastq.gz')).sample))

# Master rule for controlling workflow. 
rule all:
	input:
		"data/shared/metagenome.norm.shared",
		expand("data/shared/{group}.metagenome.norm.shared",
			group = config["metagenomeGroups"]),
		"data/catbat/metagenome.taxonomy.txt",
		expand("data/multiqc/{report}_multiqc_report.html",
			report = config["metagenomeReports"])





##################################################################
#
# Part 1: Read Processing 
#
##################################################################

# Trimming adapters and removing low quality sequences.
rule trimReadPair:
	input:
		script="code/bash/metagenomeTrimmomatic.sh",
		raw=lambda wildcards: expand("data/raw/{sample}_{readNum}_001.fastq.gz", 
            sample = wildcards.sample, readNum = config["readNum"])
	output:
		trimmedPaired1="data/trimmomatic/{sample}_R1_paired.fq.gz",
		trimmedPaired2="data/trimmomatic/{sample}_R2_paired.fq.gz",
		trimmedUnpaired="data/trimmomatic/{sample}_unpaired.fq.gz"
	conda:
		"envs/yaml/trimmomatic.yaml"
	shell:
		"bash {input.script} {input.raw}"


# Creating indices for i) mouse, ii) human, and iii) mouse + human for screen and removing host contaminants.
rule buildContaminantIndices:
	input:
		script="code/bash/metagenomeBowtie2Dependencies.sh"
	output:
		mouseIndex=expand("envs/share/bowtie2/index/m_musculus_index.{extension}",
			extension = ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
		humanIndex=expand("envs/share/bowtie2/index/h_sapiens_index.{extension}",
			extension = ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
		combinedIndex=expand("envs/share/bowtie2/index/m_musculus_h_sapiens_combined_index.{extension}",
			extension = ["1.bt2l","2.bt2l","3.bt2l","4.bt2l","rev.1.bt2l","rev.2.bt2l"])
	conda:
		"envs/yaml/bowtie2.yaml"
	shell:
		"bash {input.script}"


# Removing reads that map to host contaminants.
rule decontaminateTrimmedReads:
	input:
		script="code/bash/metagenomeDecontaminateReads.sh",
		trimmed=rules.trimReadPair.output,
		combinedIndex=rules.buildContaminantIndices.output.combinedIndex
	output:
		deconPaired1="data/bowtie2/{sample}_R1_paired_decon.fq.gz",
		deconPaired2="data/bowtie2/{sample}_R2_paired_decon.fq.gz",
		deconUnpaired="data/bowtie2/{sample}_unpaired_decon.fq.gz"
	conda:
		"envs/yaml/decon_reads.yaml"
	shell:
		"bash {input.script} {input.trimmed} {input.combinedIndex[0]}"





##################################################################
#
# Part 2: Assembly 
#
##################################################################

# Assembling reads into contigs for each sample individually.
rule assembleReads:
	input:
		script="code/bash/metagenomeSPAdes.sh",
		decon=rules.decontaminateTrimmedReads.output
	output:
		contigs="data/metaspades/{sample}/contigs.fasta"
	conda:
		"envs/yaml/spades.yaml"
	shell:
		"bash {input.script} {input.decon}"


# Combining all of the assembled contigs for all of the samples into a single fasta library.
rule combineContigs:
	input:
		script="code/bash/metagenomeCatContigs.sh",
		contigs=expand(rules.assembleReads.output.contigs,
			sample = metagenomeSamples)
	output:
		library="data/metaspades/library/contig_library.fa"
	shell:
		"bash {input.script} {input.contigs}"


# Removing sequence replicates and containments.
rule dedupeContigs:
	input:
		script="code/bash/metagenomeDedupe.sh",
		library=rules.combineContigs.output.library
	output:
		deduped="data/dedupe/contig_library_deduped.fa"
	conda:
		"envs/yaml/bbmap.yaml"
	shell:
		"bash {input.script} {input.library}"


# Removing contigs less than 1500 bp to speed up downstream processes and increase specificity.
rule filterContigs:
	input:
		script="code/bash/metagenomeSeqkitFilter.sh",
		deduped=rules.dedupeContigs.output.deduped
	output:
		filtered="data/seqkit/contig_library_deduped_filtered.fa"
	params:
		length="1500"
	conda:
		"envs/yaml/seqkit.yaml"
	shell:
		"bash {input.script} {input.deduped} {params.length}"





##################################################################
#
# Part 3: Viral Contig Curation 
#
##################################################################

# Downloading VirSorter reference database.
rule getVirsorterDatabase:
	input:
		script="code/bash/metagenomeVirSorterDependencies.sh"
	output:
		readme="envs/share/virsorter/virsorter-data/VirSorter_Readme.txt"
	shell:
		"bash {input.script}"


# Using VirSorter to predict putative viral contigs.
rule predictContigsVirsorter:
	input:
		script="code/bash/metagenomeVirSorter.sh",
		fasta=rules.filterContigs.output.filtered,
		readme=rules.getVirsorterDatabase.output.readme
	output:
		predicted="data/virsorter/virsorter_predicted_contigs.txt"
	conda:
		"envs/yaml/virsorter.yaml"
	shell:
		"bash {input.script} {input.fasta} {input.readme}"


# Downloading VirFinder eukaryotic/prokaryotic prediction model.
rule getVirfinderModel:
	input:
		script="code/bash/metagenomeVirFinderDependencies.sh"
	output:
		model="envs/share/virfinder/VF.modEPV_k8.rda"
	shell:
		"bash {input.script}"


# Using VirFinder to predict putative viral contigs.
rule predictContigsVirfinder:
	input:
		script="code/R/metagenomeVirFinder.R",
		fasta=rules.filterContigs.output.filtered,
		model=rules.getVirfinderModel.output.model
	output:
		predicted="data/virfinder/virfinder_predicted_contigs.txt"
	params:
		fpr=config["metagenomeVirfinderFpr"],
		fdr=config["metagenomeVirfinderFdr"]
	conda:
		"envs/yaml/virfinder.yaml"
	shell:
		"Rscript {input.script} {input.fasta} {input.model} {params.fpr} {params.fdr}"


# Creating master viral contig library using output from VirSorter and VirFinder.
rule curateViralContigs:
	input:
		script="code/bash/metagenomeCurateContigs.sh",
		fasta=rules.filterContigs.output.filtered,
		virsorter=rules.predictContigsVirsorter.output.predicted,
		virfinder=rules.predictContigsVirfinder.output.predicted
	output:
		curated="data/curate/viral_contig_library.fa"
	shell:
		"bash {input.script} {input.fasta} {input.virsorter} {input.virfinder}"





##################################################################
#
# Part 4: Binning 
#
##################################################################

# Preparing viral contig library for read mapping with BWA-MEM by indexing the contigs.
rule indexContigs:
	input:
		script="code/bash/metagenomeBWAIndex.sh",
		curated=rules.curateViralContigs.output.curated
	output:
		index=expand("data/bwa/contigs/index/viral_contig_library.{extension}",
			extension = ["sa","amb","ann","pac","bwt"])
	conda:
		"envs/yaml/bwa.yaml"
	shell:
		"bash {input.script} {input.curated}"


# Mapping reads to the viral contig library to assess coverage, etc.
rule mapContigs:
	input:
		script="code/bash/metagenomeBWAMap.sh",
		decon=lambda wildcards: expand("data/bowtie2/{sample}_{readNum}_paired_decon.fq.gz", 
        	sample = wildcards.sample, readNum = config["readNum"]),
		index=rules.indexContigs.output.index
	output:
		bam="data/bwa/contigs/bam/{sample}_contig_aligned.bam"
	conda:
		"envs/yaml/bwa.yaml"
	shell:
		"bash {input.script} {input.decon} {input.index[0]}"


# Binning the contigs together based on tetranucleotide frequency, coverage, and correlated abundances.
# Creates a checkpoint for populating downstream steps based on the number of bins created.
checkpoint binContigs:
	input:
		script="code/bash/metagenomeMetaBAT.sh",
		curated=rules.curateViralContigs.output.curated,
		bam=expand(rules.mapContigs.output.bam,
			sample = metagenomeSamples)
	output:
		dir=directory("data/metabat/bins") # Setting output as directory because input files are unknown
	params:
		seed=config["seed"]
	conda:
		"envs/yaml/metabat.yaml"
	shell:
		"bash {input.script} {input.curated} {params.seed} {input.bam}"


# Defining a function that pulls the names of all the contig bins from binContigs after the checkpoint finishes.
def contigBinNames(wildcards):
    checkpoint_output = checkpoints.binContigs.get(**wildcards).output.dir
    return expand("data/metabat/bins/bin.{bin}.fa",
    	bin=glob_wildcards(os.path.join(checkpoint_output, "bin.{bin}.fa")).bin)


# Downloading the CheckM references for assessing bin contamination.
rule getMarkerSets:
	input:
		script="code/bash/metagenomeCheckMDependencies.sh"
	output:
		markers="envs/share/checkM/taxon_marker_sets.tsv"
	conda:
		"envs/yaml/checkm.yaml"
	shell:
		"bash {input.script}"


# Determine classification of bins based on taxonomic tree and list of marker genes.
rule detectBinContamination:
	input:
		script="code/bash/metagenomeCheckM.sh",
		markers=rules.getMarkerSets.output.markers,
		bins=contigBinNames
	output:
		metrics="data/checkm/checkm_metrics.txt"
	conda:
		"envs/yaml/checkm.yaml"
	shell:
		"bash {input.script} {input.bins[0]}"


# Remove any bins labelled as being anything other than "root" (i.e. the bins were bacterial, etc.)
rule curateBinLibrary:
	input:
		script="code/bash/metagenomeCurateBins.sh",
		metrics=rules.detectBinContamination.output.metrics,
		bins=contigBinNames
	output:
		decon="data/metabat/library/bin_library_decon.fa",
		names="data/metabat/library/bin_library_decon_names.txt"
	shell:
		"bash {input.script} {input.metrics} {input.bins[0]}"





##################################################################
#
# Part 5: Generate Shared Bin File
#
##################################################################

# Creating indices of the bins for mapping reads against.
rule indexBins:
	input:
		script="code/bash/metagenomeBWAIndex.sh",
		decon=rules.curateBinLibrary.output.decon
	output:
		index=expand("data/bwa/bins/index/bin_library_decon.{extension}",
			extension = ["sa","amb","ann","pac","bwt"])
	conda:
		"envs/yaml/bwa.yaml"
	shell:
		"bash {input.script} {input.decon}"


# Mapping reads against bins to calculate bin abundances across samples.
rule mapBins:
	input:
		script="code/bash/metagenomeBWAMap.sh",
		deconPaired=lambda wildcards: expand("data/bowtie2/{sample}_{readNum}_paired_decon.fq.gz", 
        	sample = wildcards.sample, readNum = config["readNum"]),
		index=rules.indexBins.output.index
	output:
		bam="data/bwa/bins/bam/{sample}_bin_aligned.bam"
	conda:
		"envs/yaml/bwa.yaml"
	shell:
		"bash {input.script} {input.deconPaired} {input.index[0]}"


# Indexing mapping output bam files to prep for downstream steps.
rule indexBinBam:
	input:
		script="code/bash/metagenomeSamtoolsIndex.sh",
		bam=rules.mapBins.output.bam
	output:
		bai="data/bwa/bins/bam/{sample}_bin_aligned.bam.bai"
	conda:
		"envs/yaml/samtools.yaml"
	shell:
		"bash {input.script} {input.bam}"


# Collating the read mapping information into a table format.
rule calcBinRecruitment:
	input:
		script="code/bash/metagenomeSamtoolsIdxstats.sh",
		bam=rules.mapBins.output.bam,
		bai=rules.indexBinBam.output.bai
	output:
		idxstats="data/idxstats/{sample}_bin_aligned.idxstats.txt"
	conda:
		"envs/yaml/samtools.yaml"
	shell:
		"bash {input.script} {input.bam}"


# Combining all of the read recruitment statistics into a tidy format.
rule combineBinStats:
	input:
		script="code/R/metagenomeCombineIdxstats.R",
		idxstats=expand(rules.calcBinRecruitment.output.idxstats,
			sample = metagenomeSamples)
	output:
		combined="data/idxstats/combined_idxstats.tsv"
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.idxstats}"


# Using read recruitment data to create a shared file for all of the metagenomic bins.
rule makeBinShared:
	input:
		script="code/R/metagenomeMakeBinShared.R",
		idxstats=rules.combineBinStats.output.combined,
	output:
		shared="data/shared/metagenome.raw.shared"
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.idxstats}"


# Separating samples and controls into individual shared files.
rule splitBinShared:
	input:
		script="code/R/metagenomeSplitShared.R",
		shared=rules.makeBinShared.output.shared
	output:
		sample="data/shared/sample.metagenome.raw.shared",
		control="data/shared/control.metagenome.raw.shared"
	params:
		controls=config["metagenomeControls"]
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.shared} {params.controls}"


# Separating samples and controls into individual shared files.
rule normalizeBinShared:
	input:
		script="code/R/metagenomeNormalizeShared.R",
		shared="data/shared/{group}.raw.shared",
		idxstats=rules.combineBinStats.output.combined
	output:
		norm="data/shared/{group}.norm.shared",
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.shared} {input.idxstats}"


# Subsampling the raw metagenomic shared file to the desired number of reads before normalizing.
rule subsampleBinShared:
	input:
		script="code/R/metagenomeSubsampleShared.R",
		functions="code/R/metagenomeFunctions.R",
		shared="data/shared/{group}.metagenome.raw.shared",
		idxstats=rules.combineBinStats.output.combined
	output:
		subsample="data/shared/{group}.metagenome.subsample.norm.shared"
	params:
		subthresh=config["metagenomeSubthresh"]
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.shared} {input.idxstats} {params.subthresh}"





##################################################################
#
# Part 6: Taxonomic Classification
#
##################################################################

# Download and prepare databases needed for classifying bins with CAT/BAT.
rule getClassificationDatabases:
	input:
		script="code/bash/metagenomeCATBATDependencies.sh"
	output:
		dmnd="envs/share/catbat/database/refseq_viral.dmnd",
		nodes="envs/share/catbat/taxonomy/nodes.dmp"
	conda:
		"envs/yaml/catbat.yaml"
	shell:
		"bash {input.script}"


# Classifying bins using the RefSeq viral reference database and CAT/BAT.
rule classifyBins:
	input:
		script="code/bash/metagenomeCATBAT.sh",
		bins=contigBinNames,
		dmnd=rules.getClassificationDatabases.output.dmnd,
		nodes=rules.getClassificationDatabases.output.nodes
	output:
		metrics="data/catbat/metagenome.taxonomy.txt"
	conda:
		"envs/yaml/catbat.yaml"
	shell:
		"bash {input.script} {input.bins[0]} {input.dmnd} {input.nodes}"





##################################################################
#
# Part 7: Diversity Metrics
#
##################################################################

# Calculating diversity within samples.
rule calcBinAlphaDiversity:
	input:
		script="code/R/metagenomeAlpha.R",
		functions="code/R/metagenomeFunctions.R",
		shared=rules.splitBinShared.output.sample,
		idxstats=rules.combineBinStats.output.combined
	output:
		alpha="data/diversity/sample.metagenome.alpha.tsv"
	params:
		subthresh=config["metagenomeSubthresh"],
		alpha=config["metagenomeAlpha"]
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.shared} {input.idxstats} {params.subthresh} {params.alpha}"


# Calculating diversity between samples.
rule calcBinBetaDiversity:
	input:
		script="code/R/metagenomeBeta.R",
		functions="code/R/metagenomeFunctions.R",
		shared=rules.splitBinShared.output.sample,
		idxstats=rules.combineBinStats.output.combined
	output:
		beta="data/diversity/sample.metagenome.beta.{beta,[a-z]+}.tsv"
	params:
		subthresh=config["metagenomeSubthresh"],
		beta="{beta}"
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.shared} {input.idxstats} {params.subthresh} {params.beta}"


# Calculating PCoA ordination using beta diversity distance matrix.
rule calcBinPCoA:
	input:
		script="code/R/metagenomePCoA.R",
		beta=rules.calcBinBetaDiversity.output.beta
	output:
		stats="data/diversity/sample.metagenome.beta.{beta,[a-z]+}.pcoa.stats.tsv",
		axes="data/diversity/sample.metagenome.beta.{beta,[a-z]+}.pcoa.axes.tsv"
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.beta}"


# Calculating NMDS ordination using beta diversity distance matrix.
rule calcBinNMDS:
	input:
		script="code/R/metagenomeNMDS.R",
		beta=rules.calcBinBetaDiversity.output.beta
	output:
		stats="data/diversity/sample.metagenome.beta.{beta,[a-z]+}.nmds.stats.tsv",
		axes="data/diversity/sample.metagenome.beta.{beta,[a-z]+}.nmds.axes.tsv"
	conda:
		"envs/yaml/r.yaml"
	shell:
		"Rscript {input.script} {input.beta}"





##################################################################
#
# Part 8: Quality Control 
#
##################################################################

# Configuring Fastq Screen for checking read contamination.
rule configureContaminantReadMap:
	input:
		script="code/bash/metagenomeFastqScreenDependencies.sh",
		mouseIndex=rules.buildContaminantIndices.output.mouseIndex,
		humanIndex=rules.buildContaminantIndices.output.humanIndex
	output:
		config="envs/share/fastq_screen/fastq_screen.conf"
	conda:
		"envs/yaml/fastq_screen.yaml"
	shell:
		"bash {input.script} {input.mouseIndex[0]} {input.humanIndex[0]}"


# Checking for host (moust/human) contamination in raw reads.
rule rawContaminantReadMap:
	input:
		script="code/bash/metagenomeFastqScreen.sh",
		config=rules.configureContaminantReadMap.output.config,
		raw=lambda wildcards: expand("data/raw/{sample}_{readNum}_001.fastq.gz", 
            sample = wildcards.sample, readNum = config["readNum"])
	output:
		txt1="data/fastq_screen/raw/{sample}_R1_001_screen.txt",
		txt2="data/fastq_screen/raw/{sample}_R2_001_screen.txt"
	conda:
		"envs/yaml/fastq_screen.yaml"
	shell:
		"bash {input.script} {input.config} {input.raw}"


# Calculating read quality, adapter content, etc. for raw reads.
rule rawReadStats:
	input:
		script="code/bash/metagenomeFastQC.sh",
		raw=lambda wildcards: expand("data/raw/{sample}_{readNum}_001.fastq.gz", 
            sample = wildcards.sample, readNum = config["readNum"])
	output:
		zip1="data/fastqc/raw/{sample}_R1_001_fastqc.zip",
		zip2="data/fastqc/raw/{sample}_R2_001_fastqc.zip"
	conda:
		"envs/yaml/fastqc.yaml"
	shell:
		"bash {input.script} {input.raw}"


# Combining all raw read QC data into a single report.
rule combineRawReports:
	input:
		script="code/bash/metagenomeMultiQC.sh",
		rawMap=expand(rules.rawContaminantReadMap.output,
			sample = metagenomeSamples),
		rawStats=expand(rules.rawReadStats.output,
			sample = metagenomeSamples)
	output:
		report="data/multiqc/raw_multiqc_report.html"
	conda:
		"envs/yaml/multiqc.yaml"
	shell:
		"bash {input.script} {input.rawMap} {input.rawStats}"


# Calculating read quality, adapter content, etc. for trimmed/quality filtered reads.
rule trimmedReadStats:
	input:
		script="code/bash/metagenomeFastQC.sh",
		trimmed=lambda wildcards: expand("data/trimmomatic/{sample}_{readSet}.fq.gz", 
            sample = wildcards.sample, readSet = config["metagenomeReadSet"])
	output:
		zip1="data/fastqc/after_trimmomatic/{sample}_R1_paired_fastqc.zip",
		zip2="data/fastqc/after_trimmomatic/{sample}_R2_paired_fastqc.zip"
	conda:
		"envs/yaml/fastqc.yaml"
	shell:
		"bash {input.script} {input.trimmed}"

# Combining all QC data for trimmed/quality filtered reads into a single report.
rule combineTrimmedReports:
	input:
		script="code/bash/metagenomeMultiQC.sh",
		trimmedStats=expand(rules.trimmedReadStats.output,
			sample = metagenomeSamples)
	output:
		report="data/multiqc/trimmed_multiqc_report.html"
	conda:
		"envs/yaml/multiqc.yaml"
	shell:
		"bash {input.script} {input.trimmedStats}"


# Checking for host (moust/human) contamination in trimmed/quality filtered and decontaminated reads.
rule deconContaminantReadMap:
	input:
		script="code/bash/metagenomeFastqScreen.sh",
		config=rules.configureContaminantReadMap.output.config,
		decon=lambda wildcards: expand("data/bowtie2/{sample}_{readSet}_decon.fq.gz", 
            sample = wildcards.sample, readSet = config["metagenomeReadSet"])
	output:
		txt1="data/fastq_screen/after_decon/{sample}_R1_paired_decon_screen.txt",
		txt2="data/fastq_screen/after_decon/{sample}_R2_paired_decon_screen.txt"
	conda:
		"envs/yaml/fastq_screen.yaml"
	shell:
		"bash {input.script} {input.config} {input.decon}"


# Calculating read quality, adapter content, etc. for trimmed/quality filtered and decontaminated reads.
rule deconReadStats:
	input:
		script="code/bash/metagenomeFastQC.sh",
		decon=lambda wildcards: expand("data/bowtie2/{sample}_{readSet}_decon.fq.gz", 
            sample = wildcards.sample, readSet = config["metagenomeReadSet"])
	output:
		zip1="data/fastqc/after_decon/{sample}_R1_paired_decon_fastqc.zip",
		zip2="data/fastqc/after_decon/{sample}_R2_paired_decon_fastqc.zip"
	conda:
		"envs/yaml/fastqc.yaml"
	shell:
		"bash {input.script} {input.decon}"


# Combining all QC data for trimmed/quality filtered and decontaminated reads.
rule combineDeconReports:
	input:
		script="code/bash/metagenomeMultiQC.sh",
		deconMap=expand(rules.deconContaminantReadMap.output,
			sample = metagenomeSamples),
		deconStats=expand(rules.deconReadStats.output,
			sample = metagenomeSamples)
	output:
		report="data/multiqc/decon_multiqc_report.html"
	conda:
		"envs/yaml/multiqc.yaml"
	shell:
		"bash {input.script} {input.deconMap} {input.deconStats}"


# Calculating contig sizes, N50, etc. for all sample assemblies.
rule assemblyStats:
	input:
		script="code/bash/metagenomeQuast.sh",
		contigs=expand(rules.assembleReads.output.contigs,
			sample = metagenomeSamples)
	output:
		tsv="data/quast/report.tsv"
	conda:
		"envs/yaml/quast.yaml"
	shell:
		"bash {input.script} {input.contigs}"


# Combining assembly QC statistics into a single report.
rule combineAssemblyReports:
	input:
		script="code/bash/metagenomeMultiQC.sh",
		assemblyStats=rules.assemblyStats.output.tsv
	output:
		report="data/multiqc/assembly_multiqc_report.html"
	conda:
		"envs/yaml/multiqc.yaml"
	shell:
		"bash {input.script} {input.assemblyStats}"





##################################################################
#
# Part 9: Cleaning Up 
#
##################################################################

# Resets directory by deleting all files created by this workflow.
rule clean:
	shell:
		"""
		echo PROGRESS: Removing all workflow output.
		rm -rf envs/share/*
		rm -rf $(find data/ -mindepth 1 -maxdepth 1 -type d | grep -v "raw")
		"""
