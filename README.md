# QuaCK 
# **Qua**lity **C*ontrol **K**iosk

## V0.2.0 Release Notes
------------------------
What was updated in V0.2.0?  
* updated name from GWAS_QC_Pipeline to QuaCK
* more granularity in sex check bins reported. There are now 3 categories:  gender mixmatch, gender amiguous, and gender missing in manifest  
* rename call rate statistics in detailed report and internal report to be based on call rate instead of missing call rates  
* fixed bug with duplicate concordance name to now match on regex instead of underscore splits  
* new file output called samples_with_warnings.txt to delinieate samples that fail call rate (samples_failing_callrate_QC_details.txt) versus those with sex check issuses (samples_with_warnings.txt).  The warnings are not removed from the cleaned data set  
* more thorough clean up of temporary files in project share  
* code has been updated to be more mindful to be a little less error-prone in terms of reading in arugments such as joining paths rather than concatenating paths  
* bug fixes in code syntax  
* TO DO: fix internal report boxplots  
* TO DO: pipeline now handles strand flips and triallelic call issues for concordance  
* TO DO: SampleTable.txt is now automatically updated with call rate after removing noCall SNPs and noCall SNPs are reported  
* TO DO: internal report also gives specifics on chip failures due to poor call rate threholds  
* TO DO: duplicate concordance is now printed as a text file for every pairwise duplicate checked in addition to the average on the summary PDF  
* TO DO: better organization of a deliverable product  
* TO DO: argument added to now exclude a list of snps before proceeding with full analysis (updates new callrate in SampleTable.txt)

## Table of Contents
--------------------
1. [Overview and Purpose](#overview-and-purpose)
2. [Genome Studio Work Flow](#genome-studio-work-flow)
3. [Software Requirements](#software-requirements)
4. [Input File Requirements](#user-generated-file-requirements)
5. [Installation and Configuration](#installation-and-configuration)
6. [Running the Pipeline](#running-the-pipeline)
7. [Output and Deliverables](#output-and-deliverables)
8. [Quick Start](#quick-start-instructions)

## Overview and Purpose
------------------------
The purpose of this pipeline is to automate the first stage of genotyping quality control and create deliverables for users.  These deliverables are PDF reports regarding QC statistics on the sample and SNP level of a genotyping experiment as well as cleaned PLINK files that can be used for the round 2 of QC analysis as well as downstream analysis.  Although only a single command is used to run the pipeline once it has been installed and configured, it does have two components that run on the backend.
1. Building and Configuration of Project  
2. Running QC Pipeline 

<p align="center">
<img src="https://github.com/tbrunetti/GWAS_QC_pipeline_ver_1.0/blob/master/pre-QC-initialization-workflow.png" />
</p>
<p align="center">
<img src="https://github.com/tbrunetti/GWAS_QC_pipeline_ver_1.0/blob/master/QC-pipeline-workflow.png" />
</p>  

## Genome Studio Work Flow
--------------------------
**GS Prep Instructions and Considerations for MEGA and MEGA Custom Genotyping Chips**  
1.  Make sure all samples are included in the sample sheet (CSV)  
    * Did you include previous pilot samples that were not re-run on the full project  
    * Did you include all original and re-run samples? (due to genotyping call rate being too low initially)
2.  Sample_ID column in the sample sheet is required to be in the following format: WG[plateNumber]-DNA_[wellPosition A01-H12]_[sampleName]  --no whitespaces please!
    * Ex: WG1-DNA_B02_343523 would mean that sample ID 343523 was run on plate 1 in well position B02  
3.  When loading idats into GS, make sure none of the idats are gzipped or GS will not recognize them  
4.  Do not include samples where the well is empty or the sample failed due to wet-lab issues (ex: chip was loaded incorrectly, iScan problems, etc...)  

**Exporting Files out of GS for MEGA and MEGA Custom Genotyping Chips**  
1.  If CNV analysis is requested, perform all CNV calculation in GS using the CNV anlaysis module/plug-in  
2.  Export sample table (all columns and all samples)  
3.  Export SNP table (all columns and all snps)  
4.  Export PLINK file  
5.  Export final report (required columns are: SNP Name, Sample ID, Allele1-Top, Allele2-Top, GC Score, Sample Name, Chr, Position, Theta, R, X, Y, B Allele Freq, Log R Ratio ; if CNV is requested please make sure the following columns are added in addtion to the required columns: CNV Value, CNV Confidence)  If you are not sure or don't feel comfortable exporting specific columns, just export everything.  Additional columns will not affect the analysis.)  




## Software Requirements
------------------------
* Python version 2.7 (https://www.python.org/)
* chunkypipes (http://chunky-pipes.readthedocs.io/en/stable/getting_started.html)
* PLINK version 1.9 (https://www.cog-genomics.org/plink2)   
__*--Software Requirements that can be installed automatically--*__  
These requirements are all Python libraries/packages.  Upon installation of the pipeline, the user is given the option for the pipeline to install these dependencies automatically, assuming the user has pip installed on their system (option only available to MAC OSX and Linux operating systems, windows users need to install these libraries manually or many of these packages have EasyInstall functionality or can be easily installed from source using python setup.py install).  
  * SciPy stack, in particular the following Python packages: (https://scipy.org/)
    * pandas
    * numpy
    * matplotlib
  * Statistics (https://pypi.python.org/pypi/statistics)
  * Seaborn (http://seaborn.pydata.org/installing.html)
  * Pillow (https://pypi.python.org/pypi/Pillow/3.4.2)
  * pyFPDF (http://pyfpdf.readthedocs.io/en/latest/index.html)
  * pyPDF2 (https://pypi.python.org/pypi/PyPDF2/1.26.0)

## User Generated File Requirements
------------------------------------
There are a total of three files that the user must provide to the pipeline.
1.  Sample Table
2.  SNP Table
3.  PLINK file (PED with MAP, or BED with BIM & FAM)  

All three files can easily be generated using Illumina's free GenomeStudio software (https://www.illumina.com/techniques/microarrays/array-data-analysis-experimental-design/genomestudio.html).  However, GenomeStudio only has Windows support, therefore, unless your OS is Windows, one can create a Windows Virtual machine and run GenomeStudio there and export the results. The other option is to manually create these files using other software or by generating a tab-delimited Sample and SNP table.  The required tab-delimited headers for the sample table need to be the following (in any order):
* Sample ID
* Call Rate
* p10 GC   

The required tab-delimited headers and file for the SNP table need to be the following (in any order):
* Chr
* Name
* Cluster Sep
* AA T Mean
* AA T Dev
* BB T Mean
* BB T Dev
* AA R Mean
* AB R Mean
* BB T Mean   

Besides the Chr and Name header in the SNP file, the other file headers can have names and identifiers before and after header as long as the rest of the header name is embedded in it. (i.e. myProject_Cluster Sep would be an appropriate header because Cluster Sep is embedded in the name)

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
**Install and Configure GWAS QC Pipeline**  
After chunkypipes has been installed and initialized, the QC pipeline can be installed. First clone this repository into a directory by executing the following commands:
```
mkdir ~/my_project
cd ~/my_project
git clone https://github.com/tbrunetti/GWAS_QC_pipeline.git
cd GWAS_QC_pipeline
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
If the pipeline was installed successfully the following message will appear:
```
Pipeline run_GWAS_QC_filtering_pipeline.py sucessfully installed.
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
If yes, then the pipeline will install all the listed Python packages, if no, it is assumed the user has manually installed the listed packages.  Be aware that if the packages are not installed or improperly installed, the pipeline will not work.  **NOTE TO WINDOWS USERS!** this option can only be used if pip is installed on a Mac or Linux operating system.  Windows users must install all these manually or use a Linux virtual machine to run the pipeline.  After this, the pipeline has been successfully configured and installed.  This process will only ever need to be performed once.

## Running the Pipeline
------------------------
**Most Basic Run**  
As mentioned above there are only three files that are required to run this pipeline.  No arguments besides the arguments specifying the three files need to be called by the user.  If no arguments are set by the user except for the files, everything will be run on the default parameters and threholds.  An example of this most basic usage is shown below:
```
chunky run run_GWAS_QC_filtering_pipeline.py -sampleTable /path/to/genomeStudio_Sample_table.txt -snpTable /path/to/genomeStudio_SNPtable.txt -inputPLINK /path/to/plink/file.ped
```
The PLINK file can be generated easily using Illumina's GenomeStudio.  GenomeStudio only needs to be used to extract table information and generate a ped file.  No calculations need to be generated from GenomeStudio except for the calculations against the cluster file for SNP cluster membership and for p10GC quality for overall genotyping call rate per sample.  Additionally as long as the sample table and snp table contain the minimal headings and information listed in the file specifications, GenomeStudio is not a dependency.

___***OPTIONAL ARGUMENTS***___  
There are many options for the user to customize the pipeline for their own personal preferences.  However if no parameters and thresholds are set by the user, all arguments will be set to a default value.  These defaults are hard cut-offs recommended by Illumina or set by our lab based on what we have seen works best.

| Pipeline Argument | Default Value | Description |
| --- | --- | --- |
| --outDir | current working directory | directory to output results |
| --projectName | date time of started run | name of project |
| --arrayType | Illumina MEGA | name of genotyping chip used |
| --callrate | 0.97 | minimum sample genotying call rate for passing qc|
| --snp_callrate | 0.97 | minimum snp genotyping call rate for passing qc (Illumina hard-cut off) |
| --clusterSep | 0.30 | minimum cluster separation for passing qc (Illumina hard-cut off) |
| --AATmean | 0.30 | maximum mean of normalized AA theta values across snps (Illumina hard-cut off) |
| --AATdev | 0.06 | maximum allowable AA theta deviation threshold for passing qc (Illumina hard-cut off) |
| --BBTmean | 0.70 | minimum allowable BB theta mean threshold for passing qc (Illumina hard-cut off) |
| --BBTdev | 0.06 | maximum allowable BB theta deviation threshold for passing qc (Illumina hard-cut off) |
| --AARmean | 0.20 | minimum allowable AA mean normalized intensity for passing qc (Illumina hard-cut off) |
| --ABRmean | 0.20 | minimum allowable AB mean normalized intensity for passing qc (Illumina hard cut off) |
| --BBRmean | 0.20 | minimum allowable BB mean normalized intensity for passing qc (Illumin hard-cut off) |
| --genome_build | b37-hg19 | human genome build to use for analysis |
| --maxFemale | 0.20 | maximum estimated F coefficient for females for sex imputation |
| --minMale | 0.80 | minimum estimated F coefficient for males for sex imputation |
| --chipFailure | 1 | maximum number of sex discrepencies and failed sample missingness threshold  on a chip before considered failing |


## Output and Deliverables
--------------------------
There are 3 types of files that are output from this pipeline once it has been successfully completed: PLINK files, PDF reports, text files.  

* PLINK files - the original uncleaned PLINK files in bed format and the cleaned PLINK files in bed format that can be used for further downstream QC and analysis.
* PDF files - there are 4 PDF reports generated.  These reports aim to give the user statistics about the initial round of QC as their data progresses through the pipeline.
  * project-name_final_summary_report.pdf
  * project-name_final_detailed_report.pdf
  * project-name_final_glossary_report.pdf
  * project-name_final_internal_report.pdf
* text files - there are three text files generated.
  * snps_failing_QC_details.txt
  * samples_failing_QC_details.txt
  * md5_check_sum.txt

## Quick Start Instructions
----------------------------
Coming Soon!
