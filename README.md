# Computational Pipeline for Predicting Amino Acid Conservation Across Multiple Species

## Table of Contents
* [Aim](#aim)
* [Installation](#installation)
* [Inputs and Application](#inputs-and-application)

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
A simple and efficient way of accessing Python is by downloading the** [Anaconda Platform](https://www.anaconda.com/products/distribution) **which contains multiple relevant IDEs and allows effortless installation of any necessary Python modules and libraries.

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
-	User-prepared file *targets.fasta* with **FASTA sequences of target proteins** containing the sites of interest.
-	User-prepared file *proteomes.fasta* with **complete reference proteomes of target species** (including the origin species of the target proteins) in FASTA format across which the conservation of target sites would be assessed.
-	**MUSCLE executable file (*muscle.exe*).**
-	**Dictionary file** *Mapped_Uniprot_Species_Names.tsv* containing UniProt codes for all available species as well as their common and scientific names.
-	User-prepared CSV file *sites.csv* with **sites of interest** which must contain UniProt accession numbers of target protein sequences in the first column and positions of the sites of interest in the second column (see layout below).

The user has to then access the configuration file (*configurations.ini*) using any appropriate text editor and specify the following parameters for the conservation pipeline:
- **Origin species of the target protein sequences** (*species_of_targets* parameter). The species name must be entered using a relevant UniProt species code (if the species code is unknown, the user would refer to the UniProt database or search the pre-processed dictionary input *Mapped_Uniprot_Species_Names.tsv*).
- **Amino acid identity of target sites** (*target_amino_acid* parameter) using single-letter amino acid code. The pipeline calculates the conservation of a target amino acid at each of its positions within every query protein sequence.
- **A most likely substitution of a target amino acid** (*sub* parameter) which may not influence the site function, and which is therefore included into conservation calculation if found instead of a target amino acid within aligned sequences. If no substitution is available, the user must enter any amino acid other than the target amino acid and ignore any columns in the output referring to the substitution.
- **E-value of the BLASTp search** (*eval_thres parameter*). Any resulting BLAST hits equal to or less than the specified E-value threshold are accepted by the pipeline.

Once the parameters are specified and the configuration file is saved, the user can then run the conservation pipeline via command line by calling the main script, making sure that the location of the working directory containing all the necessary inputs is specified. 
- Location can be specified using cd command in command line, for example, *cd D:\Work\my_working_directory*
- Once the working directory is specified, Python script can be called using python command, for example, 
*python conservation_code_run_with_config.py*

It is also possible to run the script directly without using the configuration file. In this case, the user can run another version of the script *conservation_code_run_without_config.py* and run it via a relevant Python IDE.
