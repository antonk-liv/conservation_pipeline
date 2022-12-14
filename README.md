# Computational Pipeline for Predicting Amino Acid Conservation Across Multiple Species

## Table of Contents
* [Aim](#aim)
* [Installation](#installation)
* [Inputs and Application](#inputs-and-application)
* [Outputs](#outputs)
* [Troubleshooting](#troubleshooting)
* [Examples](#examples)
* [Further Analysis](#further-analysis)
* [Contact](#contact)

## Aim
This Python pipeline allows to perform the following functions in a single step:
- Search multiple query sequences against reference proteomes of selected species.
- For each query sequence, extract a top matched protein sequence from each species.
- Generate multiple sequence alignments between each query sequence and its matches.
- Determine the conservation of a specified amino acid across the aligned sequences at each of its positions in a query sequence.
- Produce a comprehensive summary output which can be mapped to target sites of interest and used to predict their conservation patterns across the selected species.

## Installation
### **Python**
In order to run the pipeline, the user must have Python programming language installed on their system. The main script can be accessed either through a relevant IDE which supports Python, via Windows/MacOS command line or via Anaconda Prompt.
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
BLAST is used in the pipeline to find protein sequence matches to each of the searched target proteins. It is required that the user has the standalone version of BLAST, [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) installed locally on their system (available from the [NIH website](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) so that the BLAST step of the pipeline can run successfully. Once installed, BLAST+ should be accessed by the pipeline automatically. However, any issues related to BLAST can be fixed by placing all the files from the downloaded BLAST+ repository into the working directory.

## Inputs and Application
All the inputs necessary for running the pipeline including the main Python script can be accessed via a [GitHub repository](https://github.com/antonk-liv/conservation_pipeline) or from the [University of Liverpool server space](https://pgb.liv.ac.uk/~antonk/).

The following eight input files must be placed into the user???s working directory before running the Python conservation pipeline:
-	**Main Python script** with the pipeline (*conservation_code_run_with_config_ver250922.py*).
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
- Location of working directory can be specified using cd command in command line, for example, *cd D:\my_working_directory*
- Python script can be called using python command, for example, *python conservation_code_run_with_config_ver250922.py*

It is also possible to run the script directly without using the configuration file. In this case, the user can run another version of the script *conservation_code_run_without_config_ver250922.py* directly using a relevant Python IDE. However, in that case the user would have to specify origin species of target proteins, target amino acid, likely substitution of the target amino acid and BLAST E-value threshold in lines 20-23 of the script.

## Outputs
### **Main output files**
Once the pipeline is run, all outputs are summarised in a newly generated folder called "ALL_OUTPUTS" where the user can find the following results:
- **Folder *Alignments_per_target*** which contains multiple sequence alignments between each protein target and its top hits in aligned FASTA (.afa) format. The alignments can be viewed directly by using relevant software such as JalView.
- **Folder *FASTA_sequences_of_top_hits_per_target*** which, for each target protein, contains its sequence in FASTA format and the sequences of its top hits.
- **Folder *Top_BLAST_hits_per_target*** which contains top BLAST hits from each species for each protein target.
- **File *alltophits.csv*** containing top BLAST hits from each species for each protein (as a single file).
- **File *blast_result_file*** containing overall raw BLAST results from BLAST+.
- **File *Conservation_of_target_X_sites.csv*** which contains summary conservation data for the target amino acid (X) sites of interest from target proteins.
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

### **Column meanings in the main conservation outputs**
The outputs containing conservation data present conservation scores in different ways and have the following columns:
| Column header | Explanation |
| ----    |:---: |
| Protein | UniProt Accession ID of a target protein |
| Peptide sequence from alignment  | Sequence of 15 amino acids containing a target site. The site is located in the middle of the sequence if -7 and +7 amino acids are available | 
| Position in alignment  | Position of target site in a multiple sequence alignment of its encompassing protein | 
| Position in target protein  | Position of target site in its protein sequence  | 
| No. of species analysed | Total number of proteomes across which conservation is determined. This number excludes the origin species of target proteins | 
| %Conservation out of total no. of proteomes | Percentage conservation score calculated out of total number of analysed species. For example, if there were 100 species of interest and the target site is conserved in 20 of them, then % conservation is 20% |
| %Conservation out of total no. of proteomes inc. substitutions  | This score takes into account the substitutions of a target amino acid when calculating the score. For example, if the target site is Ser, any conserved Ser and its most likely substitution, Thr, would be added to the numerator of a % conservation calculation  | 
| %Conservation in aligned hits  | The conservation score is calculated out of the number of aligned top protein hits from each species (i.e. the number of sequences in alignment minus 1), rather than out of the total number of target species. For example, if out of 100 target species, the alignments were made with hits from 60 species, then the denominator in the % conservation calculation would be 60 | 
| %Conservation in aligned hits inc. substitutions  | Conservation score out of the number of found hits, also taking a substitution mutation of a target amino acid into account. This means that for the % conservation calculation, numerator would be the number of conserved amino acids at a position in an alignment plus the number of likely substitutions at that position; the denominator would be the number of sequences in the alignment minus 1 | 
| No. of species aligned  | Number of target species from which a protein was matched and aligned with the target protein  | 
| Species aligned (UP codes)  | Names of aligned species given as UniProt species codes  | 
| Species aligned (common or sci names)  |  Names of aligned species given as their common or scientific names  | 
| No. of species aligned excl. gaps  | Total number of species aligned per site but only counting the species where there is no gap at a site position  | 
| %Conservation in aligned hits excl. gaps  | Conservation score out of aligned species which do not have a gap at a site position   | 
| %Conservation in aligned hits excl. gaps inc. substitutions  | Conservation score out of aligned species which do not have a gap at a site position, but also taking substitutions of a target amino acid into account and adding them to the numerator in % conservation calculation if conserved at a site position  | 
| -1 site  | Amino acid at a minus one position from the target amino acid in protein sequence  | 
| +1 site  | Amino acid at a plus one position from the target amino acid in protein sequence  | 
| -1 site position in aln.  | Position of the -1 amino acid in alignment  | 
| +1 site position in aln.  | Position of the +1 amino acid in alignment  | 
| %Conservation of -1 site in aligned hits  | Conservation of the -1 site across aligned sequences (not the total number of target species) | 
| %Conservation of +1 site in aligned hits  | Conservation of the +1 site across aligned sequences (not the total number of target species) |

## Troubleshooting
- Any issues related to BLAST (for example, blast not found) could be fixed by placing all the files from the downloaded BLAST+ repository into the working directory.
- If in the final conservation output with target sites, the positions of target sites are not found or not matched (presented by dash (-)), this means that the target sequence might be different to the one in the reference proteome for that species, even though the IDs are the same. This happens because UniProt sometimes updates the sequences and selects different isoforms as canonical, while the ID of a canonical sequence remains unchanged.

## Examples
Examples of actual inputs and outputs, as well as sample pipeline runs can be found [HERE](https://pgb.liv.ac.uk/~antonk/)

## Further Analysis
### **Analysing proteins which are not in the reference proteome of the origin species**
The pipeline was originally designed to only analyse proteins which are in a reference proteome of their origin species. However, it is possible to analyse protein targets which are not in the reference proteome by making slight adjustments to the *proteomes.fasta* input file:
- If the user already knows that a protein target is not in the reference proteome of its origin species, then the sequence of that protein target must be added to the *proteomes.fasta* input file before running the code.
- If the user does not know whether a protein target is in the reference proteome of its origin species, the user should run the code once and analyse output file *targets_not_analysed.csv*. Any protein targets which are not in the reference proteome are written into that file. The user can then extract the IDs of those proteins and create a new *targets.fasta* file containing the sequences of those proteins. Furthermore, those sequences should be added to the *proteomes.fasta* file. The code can then be re-run to analyse proteins which are not in the reference proteome.

### **Analysing additional target amino acids without re-running the full pipeline**
The pipeline was designed to analyse one amino acid at a time. However, it is possible to analyse additional amino acids after running the main code without having to re-run the full pipeline. A separate code *conservation_code_aln_section_with_config_ver250922.py* was created to isolate the section of the main code which processes the alignments and generates conservation results. The user can run the code in a similar way as running the main pipeline, but the user has to ensure that all alignments (.afa outputs) and all FASTA files containing the sequences of target proteins and their hits are also included in the working directory. As a result, the user must have the following files in the working directory to run the alignment section of the code:
- The alignment section of the code itself *conservation_code_aln_section_with_config_ver250922.py*
- Configuration file *configurations_aln_file.ini* which contains user-defined settings which must be entered by the user beforehand (Origin species of target sequences using UniProt species code, target amino acid, likely substitution, total number of proteomes excluding the origin species).
- Linker file *link_to_config_aln.py* which connects the code with the configurations file.
- Dictionary file *Mapped_Uniprot_Species_Names.tsv* containing UniProt codes for all available species as well as their common and scientific names.
-	User-prepared CSV file *sites.csv* with sites of interest which must contain UniProt accession numbers of target protein sequences in the first column and positions of the sites of interest in the second column.

As with the main pipeline, another version of the alignment section is available *conservation_code_aln_section_without_config_ver250922.py* which can be run directly using a Python IDE without interacting with a configurations file. In that case, the user has to enter the parameters within the code itself in lines 23-26.

## Contact
If you have any questions or suggestions about the conservation pipeline, feel free to contact Anton Kalyuzhnyy at A.Kalyuzhnyy@liverpool.ac.uk
