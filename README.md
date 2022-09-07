# Computational Pipeline for Predicting Amino Acid Conservation Across Multiple Species
> The Python pipeline was optimised to perform the following functions in a single step:
- Search multiple query sequences against reference proteomes of selected species.
- For each query sequence, extract a top matched protein sequence from each species.
- Generate multiple sequence alignments between each query sequence and its matches.
- Determine the conservation of a specified amino acid across the aligned sequences at each of its positions in a query sequence.
- Produce a comprehensive summary output which can be mapped to target sites of interest and used to predict their conservation patterns across the selected species.
