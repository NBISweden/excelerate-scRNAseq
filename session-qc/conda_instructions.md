# Getting started

In this workshop you will use the lab's computers, so you'll have all the necessary software installed
and ready to run everything.  

However, if you want to use your own laptop or later want to repeat these exercises, you can follow the 
instructions below to do so. These rely on conda environments, a [self-contained directory that 
you](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) can use in 
order to reproduce all your results.

Briefly, you need to:  

1. Install Conda and download the `.yaml` file
2. Create and activate the environment
3. Deactivate the environment after running your analyses

You can read more about Conda environments and other important concepts to help you make your research 
reproducible [here](https://nbis-reproducible-research.readthedocs.io/en/latest/conda/).


**Install Conda and download the environment file**

You should start by installing Conda. We suggest installing either Miniconda, or Anaconda if storage is 
not an issue. After [installing Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), 
download [environment.yaml](https://raw.githubusercontent.com/NBISweden/excelerate-scRNAseq/master/environment.yaml) 
and put it in your working folder. 

**Create and activate the environment**

In terminal, `cd` to your working folder and create an environment called `project_scRNA` from the 
`environment.yaml` file. This may take a few minutes:

```
conda env create -p project_scRNA -f environment.yaml
```
  
```
##Collecting package metadata: done
##Solving environment: done
##
##Downloading and Extracting Packages
##libcblas-3.8.0       | 6 KB      | ############################################################################# | 100%
##liblapack-3.8.0      | 6 KB      | ############################################################################# | 100%
##...
##Preparing transaction: done
##Verifying transaction: done
##Executing transaction: done
###
### To activate this environment, use
###
###     $ conda activate /my/file/path/project_scRNA
###
### To deactivate an active environment, use
###
###     $ conda deactivate
```

Then activate it:
```
conda activate ./project_scRNA
```

From this point on you can run any of the contents from the course. For instance, you can directly launch RStudio by 
typing `rstudio`.

**Deactivate the environment**
After you've ran all your analyses, deactivate the environment by typing `conda deactivate`.
