# Computational Pipeline for Predicting Amino Acid Conservation Across Multiple Species

## Table of Contents
* [Aim](#aim)
* [Installation](#installation)

## Aim
**This Python pipeline allows to perform the following functions in a single step:**
- Search multiple query sequences against reference proteomes of selected species.
- For each query sequence, extract a top matched protein sequence from each species.
- Generate multiple sequence alignments between each query sequence and its matches.
- Determine the conservation of a specified amino acid across the aligned sequences at each of its positions in a query sequence.
- Produce a comprehensive summary output which can be mapped to target sites of interest and used to predict their conservation patterns across the selected species.

## Installation
#### **Python**
**In order to run the pipeline, the user must have Python programming language installed on their system. The main script can be accessed either through a relevant IDE which supports Python, via Windows/MasOS command line or via Anaconda Prompt.**

**A simple and efficient way of accessing Python is by downloading the** [Anaconda Platform](https://pages.github.com/) **which contains multiple relevant IDEs and allows effortless installation of any necessary Python modules and libraries.**

###### **Python Modules and Libraries**
**Once Python is installed, the user must ensure that all the modules and libraries used in the pipeline are installed and updated.**

| Python module/library | Version used | Installation and more info |
| ------------- |:-------------:| -----------------------------------------------------------:|
| Numpy         | 1.21.5        | [Numpy more info](https://numpy.org/)                       |
| BioPython     | 1.74          | [BioPython](https://biopython.org/)                         |
| CSV           | 1.0           | [CSV more info](https://docs.python.org/3/library/csv.html) |
| Re            | 2.2.1         | [Re more info](https://docs.python.org/3/library/re.html)   |
| Pandas        | 1.3.5         | [Pandas more info](https://pandas.pydata.org/)              |

**If [Anaconda](https://pages.github.com/) is available on the system, the user can install all libraries via Anaconda Navigator under *Environments* or via Anaconda Prompt using a relevant installation command for a given module/library.**

###### **BLAST+**
**BLAST is used in the pipeline to find protein sequence mathces to each of the searched target proteins. It is required that the user has the standalone version of BLAST, [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) installed locally on their system (available from the [NIH website](https://blast.ncbi.nlm.nih.gov/Blast.cgi) so that the BLAST step of the pipeline can run successfully. Once installed, BLAST+ should be accessed by the pipeline automatically. However, any issues related to BLAST can be fixed by placing all the files from the downloaded BLAST+ repository into the working directory.**
