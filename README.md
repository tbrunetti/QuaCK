# GWAS_QC_pipeline_ver_1.0

## Overview and Purpose
------------------------
The purpose of this pipeline is to automate the first stage of genotyping quality control and create deliverables for users.  These deliverables are PDF reports regarding QC statistics on the sample and SNP level of a genotyping experiment as well as cleaned PLINK files that can be used for the round 2 of QC analysis as well as downstream analysis.  Although only a single command is used to run the pipeline once it has been installed and configured, it does have two components that run on the backend.
1. Building and Configuration of  Project  
2. Running QC Pipeline 

<p align="center">
<img src="https://github.com/tbrunetti/GWAS_QC_pipeline_ver_1.0/blob/master/pre-QC-initialization-workflow.png" />
</p>
<p align="center">
<img src="https://github.com/tbrunetti/GWAS_QC_pipeline_ver_1.0/blob/master/QC-pipeline-workflow.png" />
</p>

## Software Requirements
------------------------
* Python version 2.7.6
* chunkypipes (http://chunky-pipes.readthedocs.io/en/stable/getting_started.html)
* PLINK version 1.9 (https://www.cog-genomics.org/plink2)
__*Software Requirements that can be installed automatically*__
These requirements are all Python libraries/packages.  Upon installation of the pipeline, the user is given the option for the pipeline to install these dependencies automatically, assuming the user has pip installed on their system (option only available to MAC OSX and Linux operating systems, windows users need to install these libraries manually).
* SciPy stack, in particular the following Python packages:
  * pandas
  * numpy
  * matplotlib
* Statistics
* Seaborn
* Pillow
* pyFPDF
* pyPDF2

## User Generated File Requirements
------------------------------------

## Installation and Configuration
----------------------------------
First chunkypipes should be downloaded from https://pypi.python.org/pypi/ChunkyPipes and follow the installation instructions.  Please note that Mac OS and Linuz users can install chunkypipes via the command line by the following command:
'''
pip install chunkypipes
'''
Note that sudo privledges will be required for system-wide installs.
## Running the Pipeline
------------------------

___***OPTIONAL ARGUMENTS***___

## Ouput and Deliverables
--------------------------
