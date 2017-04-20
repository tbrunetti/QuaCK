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
__*--Software Requirements that can be installed automatically--*__  
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
**Download and Install chunkypipes**
First chunkypipes should be downloaded from https://pypi.python.org/pypi/ChunkyPipes and follow the installation instructions.  Please note that Mac OS and Linux users can install chunkypipes via the command line by the following command:  
```
pip install chunkypipes
```

Note that sudo privileges will be required for system-wide installs.  Once chunkypipes has been downloaded and installed run the following command:
```
chunky init
```
This will initialize the directory for which configured pipelines will be stored.  This only needs to be performed once for any chunkypipes usage.  If chunkypipes was successfully installed the follow message will appear:
```
> ChunkyPipes successfully initialized at /home/user
```
**Install GWAS QC Pipeline**
After chunkypipes has been installed and initialized, the QC pipeline can be installed. First clone this repository into a directory by executing the following commands:
```
mkdir ~/my_project
cd ~/my_project
git clone https://github.com/tbrunetti/GWAS_QC_pipeline_ver_1.0.git
cd GWAS_QC_pipeline_ver_1.0
```
Assuming chunkypipes has been installed correctly, run the following command to configure the pipeline:
```
chunky configure run_GWAS_QC_filtering_pipeline.py
```
This will prompt the user for the following information:  It is CRITICAL that the full file path is written out at each prompt.
```
Full path to PLINK executable (must be version >=1.9) []:
```
Upon successful configuration the following message will appear:
```
Configuration file successfully written.
```
Now the pipeline is successfully configured.  It will never have to be configured again unless the location of the executable for PLINK has moved to another location.   
The last thing to do is to install the pipeline now that the configuration file has been saved.  To install the pipeline, run the following command:
```
chunky install run_GWAS_QC_filtering_pipeline.py
```
This will ask the user the following question:
```
```
If the pipeline has been previously installed it will ask the user if they would like to overwrite the existing pipeline.  This will only ever be required if changes have been made to the code.
```
Pipeline run_GWAS_QC_filtering_pipeline.py is already installed, overwrite? [y/n] 
```
Upon successful installation the following message should appear:
```
Pipeline run_GWAS_QC_filtering_pipeline.py successfully installed.
```
Following this message an additional prompt appears on the command line:
```
Attempting to install the following dependencies:
pandas
matplotlib
fpdf
Pillow
seaborn
pypdf2

Proceed with dependency installation? [y/n] 

``` 
If yes, then the pipeline will install all the listed Python packages, if no, it is assumed the user has manually installed the listed packages.  If not, the pipeline will not work.  Note, this option can only be used if pip is installed on a Mac or Linux operating system.  Windows users must install all these manually or use a Linux virtual machine to run the pipeline.  After this, the pipeline has been successfully configured and installed.  This process will only ever need to be performed once.

## Running the Pipeline
------------------------

___***OPTIONAL ARGUMENTS***___

## Ouput and Deliverables
--------------------------
