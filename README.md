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
**In order to run the pipeline, the user must have Python programming language installed on their system. The main script can be accessed either through a relevant IDE which supports Python, via Windows/MasOS command line or via Anaconda Prompt.**

**A simple and efficient way of accessing Python is by downloading the** [Anaconda Platform](https://pages.github.com/) **which contains multiple relevant IDEs and allows effortless installation of any necessary Python modules and libraries.**

**Once Python is installed, the user must ensure that all the modules and libraries used in the pipeline are installed and updated (see the table below).**

| Python module/library | Version used | Installation and more info |
| ------------- |:-------------:| -----:|
| Numpy      | 1.21.5 | $1600 |
| BioPython      | 1.74      |   $12 |
| CSV | 1.0      |    $1 
| Re | 2.2.1      |    $1 
| Pandas | 1.3.5      |    $1 

