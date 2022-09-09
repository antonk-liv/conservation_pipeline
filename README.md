# Computational Pipeline for Predicting Amino Acid Conservation Across Multiple Species

## Table of Contents
* [Aim](#aim)
* [Installation](#installation)
* [Inputs and Application](#inputs-and-application)
* [Outputs](#outputs)

## Aim
This Python pipeline allows to perform the following functions in a single step:
- Search multiple query sequences against reference proteomes of selected species.
- For each query sequence, extract a top matched protein sequence from each species.
- Generate multiple sequence alignments between each query sequence and its matches.
- Determine the conservation of a specified amino acid across the aligned sequences at each of its positions in a query sequence.
- Produce a comprehensive summary output which can be mapped to target sites of interest and used to predict their conservation patterns across the selected species.

## Installation
### **Python**
In order to run the pipeline, the user must have Python programming language installed on their system. The main script can be accessed either through a relevant IDE which supports Python, via Windows/MasOS command line or via Anaconda Prompt.
A simple and efficient way of accessing Python is by downloading the **[Anaconda Platform](https://www.anaconda.com/products/distribution)** which contains multiple relevant IDEs and allows effortless installation of any necessary Python modules and libraries.

### **Python Modules and Libraries**
Once Python is installed, the user must ensure that all the modules and libraries used in the pipeline are installed and updated.

| Python module/library | Version used | Installation and more info |
| ------------- |:-------------:| -----------------------------------------------------------:|
| Numpy         | 1.21.5        | [Numpy more info](https://numpy.org/)                       |
| BioPython     | 1.74          | [BioPython](https://biopython.org/)                         |
| CSV           | 1.0           | [CSV more info](https://docs.python.org/3/library/csv.html) |
| Re            | 2.2.1         | [Re more info](https://docs.python.org/3/library/re.html)   |
| Pandas        | 1.3.5         | [Pandas more info](https://pandas.pydata.org/)              |

If [Anaconda](https://www.anaconda.com/products/distribution) is available on the system, the user can install all libraries via Anaconda Navigator under *Environments* or via Anaconda Prompt using a relevant installation command for a given module/library.

### **BLAST+**
BLAST is used in the pipeline to find protein sequence mathces to each of the searched target proteins. It is required that the user has the standalone version of BLAST, [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) installed locally on their system (available from the [NIH website](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) so that the BLAST step of the pipeline can run successfully. Once installed, BLAST+ should be accessed by the pipeline automatically. However, any issues related to BLAST can be fixed by placing all the files from the downloaded BLAST+ repository into the working directory.

## Inputs and Application
All the inputs necessary for running the pipeline including the main Python script can be accessed via a [GitHub repository](https://github.com/antonk-liv/conservation_pipeline)

The following eight input files must be placed into the user’s working directory before running the Python conservation pipeline:
-	**Main Python script** with the pipeline (*conservation_code_run_with_config.py*).
-	**Configuration file** with user-specified settings (*configurations.ini*).
-	**Linker Python file** (*link_to_config.py*) which connects the main conservation code with the configuration file. 
-	User-prepared file *targets.fasta* with **FASTA sequences of target proteins** containing the sites of interest. The sequences of protein targets can be searched and downloaded in FASTA format from UniProt.
-	User-prepared file *proteomes.fasta* with **complete reference proteomes of target species** (including the origin species of the target proteins) in FASTA format across which the conservation of target sites would be assessed. The proteomes can be downloaded from UniProt.
-	**MUSCLE executable file (*muscle.exe*).** The file is available in the GitHub directory but can also be downloaded from the [MUSCLE website](https://drive5.com/muscle/).
-	**Dictionary file** *Mapped_Uniprot_Species_Names.tsv* containing UniProt codes for all available species as well as their common and scientific names.
-	User-prepared CSV file *sites.csv* with **sites of interest** which must contain UniProt accession numbers of target protein sequences in the first column and positions of the sites of interest in the second column (see layout below).

| Protein | Site |
| ----    |:---: |
| A1A4V9  | 142  |
| A1A4V9  | 328  | 
| A1A4Y4  | 42   | 
| A6NIE6  | 153  | 
| O60361  | 29   | 
| O60361  | 110  |

The user has to then access the configuration file (*configurations.ini*) using any appropriate text editor and specify the following parameters for the conservation pipeline:
- **Origin species of the target protein sequences** (*species_of_targets* parameter). The species name must be entered using a relevant UniProt species code (if the species code is unknown, the user would refer to the UniProt database or search the pre-processed dictionary input *Mapped_Uniprot_Species_Names.tsv*).
- **Amino acid identity of target sites** (*target_amino_acid* parameter) using single-letter amino acid code. The pipeline calculates the conservation of a target amino acid at each of its positions within every query protein sequence.
- **A most likely substitution of a target amino acid** (*sub* parameter) which may not influence the site function, and which is therefore included into conservation calculation if found instead of a target amino acid within aligned sequences. If no substitution is available, the user must enter any amino acid other than the target amino acid and ignore any columns in the output referring to the substitution.
- **E-value of the BLASTp search** (*eval_thres parameter*). Any resulting BLAST hits equal to or less than the specified E-value threshold are accepted by the pipeline.

Once the parameters are specified and the configuration file is saved, the user can then run the conservation pipeline via a standard command prompt or Anaconda Prompt by calling the main script and making sure that the location of the working directory containing all the necessary inputs is specified. 
- Location of working directory can be specified using cd command in command line, for example, *cd D:\Work\my_working_directory*
- Python script can be called using python command, for example, *python conservation_code_run_with_config.py*

It is also possible to run the script directly without using the configuration file. In this case, the user can run another version of the script *conservation_code_run_without_config.py* directly using a relevant Python IDE. However, in that case the user would have to specify origin species of target proteins, target amino acid, likely substitution of the target amino acid and BLAST E-value threshold in lines 20-23 of the script.

## Outputs
### **Main output files**
Once the pipeline is run, all outputs are summarised in a newly generated folder called "ALL_OUTPUTS" where the user can find the following results:
- **Folder *Alignments_per_target*** which contains multiple sequence alignments between each protein target and its top hits in aligned FASTA (.afa) format. The alignments can be viewed directly by using relevant software such as JalView.
- **Folder *FASTA_sequences_of_top_hits_per_target*** which, for each target protein, contains its sequence in FASTA format and the sequences of its top hits.
- **File *Conservation_of_target_sites.csv*** which contains summary conservation data for the target amino acid sites of interest from target proteins.
- **File *X_conservation_in_TARGET_SPECIES.csv*** where X is the target amino acid and TARGET_SPECIES is the species code of the origin species from which the targets came from. This is a summary file containing conservation data for every target amino acid in each target protein.
- **File *targets_not_analysed.csv*** containing protein targets which were not analysed and the reasons for their exclusion (see below for more details).

### **Protein targets which were excluded from the analysis and the reasons for their exclusion**
Our conservation pipeline was optimised to process as many target proteins as possible. However, some targets cannot be analysed by the pipeline due to one of the following reasons:
- A target produces no hits in a BLAST search at all, likely due to its sequence being too short or poorly annotated.
- A target produces no significant hits in a BLAST search that meet the set E-value threshold, likely due its sequence being unique to its origin species.
- A target is not found in the reference proteome of its origin species, likely due to poor sequence annotation or because the target is an isoform of a canonical protein from a reference proteome. This is detected when an identical match to the query protein is not found in a BLAST search when searched against the reference proteome of its origin species.
- Multiple sequence alignments are not generated between the target sequence and its top hits. This is likely to happen when very large sequences of >30,000 base pairs are being aligned.
- A target protein does not have a target amino acid in its sequence.

Those targets are identified and summarised in the output file *targets_not_analysed.csv*.
