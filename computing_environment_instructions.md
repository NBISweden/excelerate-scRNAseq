# Using the course computing environments

You can use the classroom computers where all the necessary software has been installed as a Conda environment. In order to activate the environment
1. Log in to the computer as cscuser
2. Open the terminal (right click on desktop)
3. Go to the home directory
4. Type `source setup-single-rna-env.sh` and the prompt should change to `(single-rna-env)`
5. Launch Rstudio by typing `rstudio &`
6. All the course data is available in the folder `scrna-seq2019`
7. When the course is over you can deactivate the Conda environment by typing `conda deactivate`

If you would like to have the same software environment on your own computer after the course, please read the [Conda installation instructions](conda_instructions.md). 

The course environment including the data is also available as a virtual machine image in CSC's cPouta cloud. You can access Rstudio server on the following 10 instances, ask Eija for login details if needed:

* 193.167.189.100:8787
* 195.148.31.253:8787
* 86.50.169.110:8787
* 86.50.169.17:8787
* 86.50.169.11:8787
* 86.50.169.46:8787
* 86.50.169.107:8787
* 86.50.169.121:8787
* 86.50.169.117:8787
* 86.50.169.89:8787
