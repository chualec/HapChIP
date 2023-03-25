Software Requirement Specification
For: hichip
by: Alec

1. Introduction
- Purpose of this Document:
- HiChIP data is a new form of data that combines both HiC chromatin-capture approaches with protein pull-downs to obtain highly specific long-distance reads for specific markers including CTCF, H3K4me3, and H3K27ac. These new sequencing technologies allow us to look at chromatin interactions in a spacially-resolved way, however many tools for efficient standardization and analysis of HiChIP data is still lacking. This document hopes to give an introduction in the why looking at chromatin interactions in a haplotype-resolved way can provide greater context for chromatin interactions, how to install and use this tool, and how to interprate the results from this tool.
- Scope of the tool:
	- This tool will perform rudimentary quality control including:
	- quantifying percentage of high quality reads mapping uniquly to distal locations.
	- quantifying percent of trans-chromosome reads.
	- perform plotting of quality control
- Additionally this tool will perform basic analysis including
	- estimating differnetially bound locations between haplotypes.
	- Create rudimentary map of significant distal cis-interactions.
	- Create rudimentary map of significant trans-interactions.
  
2. Overall description
- The inputs for this tool will include:
	- fastq files or bam files of HiChIP data.
	- reference file for alignment. (can be downloaded)
	- gtf file for gene quantification (can be downloaded)
	- eQTL file for hotspot detection (optional, can be downloaded)
b. The outputs for this tool will include:
	- VCF file containing haplotype-resolved variants.
- 3 bam files containing different interactions observed in the data:
	- bam file of haplotype 1 cis-interaction
	- bam file of haplotype 2 cis-interaction
	- bam file of trans-interaction
3. Prerequisites:
- Required packages include: Conda environment, snakemake tools.


Installation
``` sh
# Clone the repository
git clone https://github.com/chualec/HapChIP

# Move to the folder
cd HapChIP

# Install necessary packages using conda
conda create env --name env --file environment.yaml
```


example run
``` sh
python3 HapChIP/HapChIP.py -b example_input/example.bam -v example_input/example.vcf.gz -o ./example_output
```


  
A unit test has been provided with a bam file and vcf file covering chr13 region.
The output will contain 5 files:
- 4 bam files containing reads split into either haplotype 1, haplotype 2, unknown haplotypes, or conflicting haplotypes
- One log file containing overall summary
