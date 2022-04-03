Welcome to the SiFoN repository (Seeker of intersecting Functions of Noncoding mutations). 

SiFoN is a framework built on top of [Sei](https://github.com/FunctionLab/sei-framework) that can be used to to understand the combined impact of interacting mutations, to visualize the impact of mutations in their genomic context, and to combine clinical information with functional predictions. SiFoN includes a pipeline of analytical visualization tools that prioritizes combinations of SNPs predicted to have causal impacts on disease phenotypes. 

The source code is available in `src/SiFoN`. 
* `tutorials` folder includes Jupyter notebook walkthroughs for each module. These walkthroughs provide examples of use cases and test data for each module. 
* `tutorials/test_input_data` includes input data for walkthroughs. 
* `tutorials/test_output_plots` includes plot outputs for the walkthroughs. 
* `model_data` includes Sei sequence class names and chromatin profile names. 
* `requirments.txt` includes all packages required to run SiFoN. 

Note that the following file structure is assumed for SiFoN:
``` bash
├── src/
│   ├── SiFoN
├── directory/
│   ├── notebook/script that imports SiFoN
```

Also note that there are three files required for the tutorial that are available on [DropBox](https://www.dropbox.com/scl/fo/tz61l2a1kmxma4p10yc7d/h?dl=0&rlkey=jwz5avk5ara6im747x2od697o). If you are running these tutorials, then please download these files and include them in your copy of `tutorials/test_input_data`. You can also use your own test data (if it is in the correct format) by including the data in `tutorials/test_input_data` and changing the notebook to read in the data.
* `output.bed` is used for `tutorials/GWAS_enrichment.ipynb` 
* `chr12_115750000_116000000_sequence_class_scores.npy` is used for `tutorials/smoothing.ipynb`
* `chr10_89580225_89633389_diffs` is used for `tutorials/chrom_viz.ipynb`

Please visit the [SiFoN documentation page](https://bmacedo-lgtm.github.io/SiFoN/) for a detailed description of each module and their functions.
