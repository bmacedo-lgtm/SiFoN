Welcome to the SiFoN (Seeking intersecting Functions of Noncoding mutations) repository. 
SiFoN is framework built on top of [Sei](https://github.com/FunctionLab/sei-framework) that can be used to understand the combined 
impact of interacting mutations, visualize the impact of mutations in their genomic context at various nucleotide scales, 
and combine clinical information with functional predictions. 
SiFoN includes a pipeline of analytical visualization tools that prioritizes combinations of SNPs predicted to have 
causal impacts on disease phenotypes. SiFoN includes a plotting library that provides efficient and user-friendly visualizations 
that highlight what regulatory roles are important to a risk loci of interest and what SNPs have the highest likelihood of causing regulatory dysfunction. 

`requirements.txt` contains a list of packages that are required to run Sei.
`test.ipynb` is a walkthrough of all the functions available in SiFoN with accompanying test data and outputs. Large outputs are available in the `test_output_plots` folder.

The SiFoN package includes four modules:
* construst_FASTA.py : functions to construct FASTA files with haplotype sequences.
* chrom_viz.py : visualizations to study the affected chromatin profiles of SNPs.
* viz.py : visualizations of sequence class predictions along genomic sequences.
* GWAS_enrichment.py : visualizations of the enrichment of sequence class annotations in GWAS risk loci.
