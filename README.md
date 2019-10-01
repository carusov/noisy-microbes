# noisy-microbes
Benchmarking of 16S clustering methods on low biomass microbiome data

This repo contains the code used for the analysis presented in [*Performance of microbiome sequence inference methods in environments of varying biomass*](https://msystems.asm.org/content/4/1/e00163-18). This study benchmarked selected 16S rRNA clustering methods on high and low microbial biomass 16S rRNA gene sequencing data. 

##Motivation
Sequence clustering is one of the initial and critical processing steps in 16S rRNA gene sequencing studies. Historically, this was completed using operational taxonomic unit (OTU) clustering methods. Recently, newer methods have been developed which overcome come limitaions of OTU clustering algorithms. These methods allow for clustering / resolution of exact sequence variants and are commonly refered to as amplicon sequence variant (ASV) methods. 


## Scripts
### Sample processing
Mock community samples are processed by a series of nested scripts, as shown here:

- **process_data.sh** - top level
	- **process_high_biomass.sh** - high biomass samples
	- **process_dilution.sh** - dilution series samples
		- **merge_and_filter.sh** - merge read pairs, quality filter
	    - **rmd2r.R** - convert .Rmd to .R script
        - **fastq_truncate.Rmd** - trim low quality tails
    - **run_all_pipelines.sh** - execute clustering/sequence inference
        - **run_5_pipelines.sh** - execute 5 of 6 pipelines
		    - **uclust_pipeline.sh** - execute UCLUST algorithm
			- **uparse_pipeline.sh** - execute UPARSE algorithm
			- **unoise_pipeline.sh** - execute UNOISE algorithm
			- **med_pipeline.sh** - execute MED algorithm
			- **deblur_pipeline.sh** - execute Deblur algorithm
        - **dada2_pipeline.Rmd** - execute DADA2 algorithm

### Analysis
Analysis of clustering/inference results is contained in two Rmd scripts, with the exception of the blast routines.

- **high_biomass_analysis.Rmd** - analyze high biomass results
- **dilution_analysis.Rmd** - analyze dilution series results
	- **blast_all.sh** - run BLAST on all fastas
		- **blast_seqs.sh** - run BLAST on a fasta file
