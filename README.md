Welcome to the SiFoN repository (Spotter of intersecting Functions of Noncoding mutations). 

SiFoN is a framework built on top of [Sei](https://github.com/FunctionLab/sei-framework) that can be used to to understand the combined impact of interacting mutations, to visualize the impact of mutations in their genomic context, and to combine clinical information with functional predictions. SiFoN includes a pipeline of analytical visualization tools that prioritizes combinations of SNPs predicted to have causal impacts on disease phenotypes. SiFoN includes a plotting library that provides efficient and user-friendly visualizations that highlight what regulatory roles are important to a risk loci of interest and what SNPs have the highest likelihood of causing regulatory dysfunction. 

The source code is available in `src/SiFoN`. 
* `tutorials` folder includes Jupyter notebook walkthroughs for each module. These walkthroughs provide examples of use cases and test data for each module. 
* `tutorials/test_input_data` includes input data for walkthroughs. 
* `tutorials/test_output_plots` includes plot outputs for the walkthroughs. 
* `model_data` includes Sei sequence class names and chromatin profile names. 
* `requirments.txt` includes all packages required to run SiFoN. 

Note that the following file structure is assumed for SiFoN:
├── src/

│   ├── SiFoN

├── directory/

│   ├── notebook/script that imports SiFoN

Please visit the [SiFoN documentation page](https://bmacedo-lgtm.github.io/SiFoN/) for a detailed description of each module and their functions.
