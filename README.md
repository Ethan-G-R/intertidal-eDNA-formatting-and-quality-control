# Formatting and quality control of eDNA metabarcoding data (intertidal dataset)

*Dina-Leigh Simons*

This GitHub repository contains R scripts and data for data formatting and performing quality control checks on an intertidal eDNA metabarcoding dataset (CO1 and 18S). This pipeline aims to bridge outputs from the bioinformatics processing (in this case, DADA2) to ecological community analyses. This repository reproduces part of the analysis pipeline for the following publications:
- Simons et al (2025), Characterising rocky intertidal biodiversity using environmental DNA metabarcoding from local to national scales.

## Sources of data used here

We use outputs from a bioinformatics pipeline for processing eDNA metabarcoding sequence data (Illumina) developed by the NERC Environmental Omics Facility (NEOF) using the University of Sheffield's computing cluster (BESSEMER), available [here](https://github.com/khmaher/HPC_dada2). We also use metadata corresponding to intertidal samples described here.

## Repository content

- Input_Data: Required data inputs from the bioinformatics pipeline described above, including BLAST taxonomic assignments (.txt), ASV information (.tsv and .fasta), and metadata (.csv). These input files correspond to multiple sampling regions in the UK.
- Processed_Data: Processed datasets generated in this pipeline.
- eDNA_Formatting_and_Quality_Control.qmd: The [Quarto](https://quarto.org/docs/get-started/hello/rstudio.html) Markdown file.
- Figures: Image files generated by eDNA_Formatting_and_Quality_Control.qmd.
- eDNA_Formatting_and_Quality_Control.RProj: The R project file. Open this project locally in RStudio.
- eDNA_Formatting_and_Quality_Control.html: The output from a successful knit of eDNA_Formatting_and_Quality_Control.qmd.

## Using this repository

The pipeline can executed by running the **`eDNA_Formatting_and_Quality_Control.qmd`** Quarto Markdown file.

#### Step 1: Check R and RStudio is installed

The code was tested against *R* 4.4.2. *R* is available for download [here](https://www.r-project.org/) and RStudio is available for download [here](https://posit.co/download/rstudio-desktop/).

#### Step 2: Download GitHub repository as an RStudio project

Download the repository via GitHub, either via RStudio or in the terminal.

#### Step 3: Work through 'eDNA_Formatting_and_Quality_Control.qmd'

Run the code chunks in the console, or knit to an .html or .pdf to see the full output.

## License

Check licenses
