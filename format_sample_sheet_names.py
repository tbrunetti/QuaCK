import pandas as pd

# spreadsheet order with header names must be spelled out exactly for genomeStudio:
# Sample_ID, Sample_Name, Sample_Plate, Sample_Well, SentrixBarcode_A, SentrixPosition_A, Gender, Sample_Group, Replicates, Parent1
# only Sample_ID, SentrixBarcode_A, and SentrixPosition_A are required columns for use in genomeStudio, Gender is required for the GWAS QC Pipeline
# an added [Header] section, [Manifest] section, and [Data] section must also be added (all the stuff from this code goes under data section)

def reformat(file):

	# Sample_ID, Sample_Plate and Sample_Well are required for QC Pipeline to run
	sample_sheet = pd.read_excel(file, sheetname="Plates_1-10", header=0)
	sample_sheet['Sample_ID'] = 'WG'+sample_sheet['Sample_Plate'].astype(str) + '-DNA_' + sample_sheet['Sample_Well'].astype(str)+ '_'+sample_sheet['Sample_Name'].astype(str)
 	sample_sheet.to_csv(file[:-5] + '_reformatted.csv', index=False, sep=',') # MUST BE COMMA SEPARATED FOR GENOMESTUDIO!

if __name__ == '__main__':
	#file = '/home/tonya/Downloads/Divers_project_master_template_5_2017.xlsx'
	file='/home/tonya/Downloads/Divers_Master_Sample_Sheet_allPlates.xlsx'
	reformat(file)


'''
Illumina sample sheet format for genomeStudio This file must be placed in same directory as SentrixBarcode_A file directory locations
All the following info needs to be added manually to the newly created csv file before importing into GenomeStudio
[Header]
Investigator Name,Name
Project Name,Project1
Experiment Name,Plate1-10
Date,5-June-16

[Manifests]
A,Multi-EthnicGlobal_A1.bpm

[Data] 
Sample_ID,SentrixBarcode_A,SentrixPosition_A,,,,,, (this part including header name is taken care of above)

'''