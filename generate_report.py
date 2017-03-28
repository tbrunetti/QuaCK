from fpdf import FPDF
import pandas
import matplotlib.pyplot as plt
import statistics as stats
import collections
import numpy as np


def thresholds_and_parameters(callrate):
	pass;

def illumina_sample_overview(inputFile, pdf, project, callrate):
	samples_to_remove_text = open('samples_to_remove.txt', 'w')
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)

	sample_qc_table = pandas.read_table(inputFile)
	total_samples = len(list(sample_qc_table['Sample ID']))
	# retrieve sample ids of those with missing call rage less than call rate parameter provided by user; default = 0.991
	samples_to_remove = list(sample_qc_table[sample_qc_table['Call Rate'] < callrate]['Sample ID'])
	basic_call_stats = [stats.median(sample_qc_table['Call Rate']), stats.mean(sample_qc_table['Call Rate']), stats.stdev(sample_qc_table['Call Rate']), min(sample_qc_table['Call Rate']), max(sample_qc_table['Call Rate'])]
	pdf.cell(0, 30, "Sample Quality", 0, 1, 'L')
	pdf.set_font('Arial', '', 16)
	pdf.cell(0, 8, "Total samples analyzed:  "+str(total_samples), 0, 1, 'L')
	pdf.cell(0, 8, "Number of samples passing missing call rate threshold:  " + str(total_samples-len(samples_to_remove)), 0, 1, 'L')	
	pdf.cell(0, 8, "Median call rate:  "+ str(basic_call_stats[0]), 0, 1, 'L')
	pdf.cell(0, 8, "Mean call rate:  "+ str(basic_call_stats[1]), 0, 1, 'L')
	pdf.cell(0, 8, "Standard deviation call rate:  "+ str(basic_call_stats[2]), 0, 1, 'L')
	pdf.cell(0, 8, "Minimum call rate:  "+ str(basic_call_stats[3]), 0, 1, 'L')
	pdf.cell(0, 8, "Maximum missing call rate:  "+ str(basic_call_stats[4]), 0, 1, 'L')
	
	# ethnicity and race breakdown
	sample_qc_table['RACE'].replace(np.nan, 'NaN', regex=True, inplace=True) # converts pandas nans into strings (by default they are ints)
	sample_qc_table['RACE'].value_counts(dropna=False).plot(kind='pie')
	plt.axis('equal')
	plt.title('Ethnicity')
	plt.savefig('ethinic_breakdown.png', bbox_inches='tight')
	pdf.image('ethinic_breakdown.png', 0, 200, 90, 70)
	
	all_samples_ethnicity = dict(sample_qc_table['RACE'].value_counts(dropna=True)) # creates the value counts of ethnicity into dictionary for easy PDF writing key=eth; value=total samples


	store_removal_ethnicity = []
	for i in samples_to_remove:
		 store_removal_ethnicity.append(list(sample_qc_table.loc[sample_qc_table['Sample ID'] == i]['RACE']))
	store_removal_ethnicity = [store_removal_ethnicity[i][0] for i in range(0, len(store_removal_ethnicity))]
	ethnicity_removals = collections.Counter(store_removal_ethnicity)
	
	pdf.set_font('Arial', 'B', 20)
	pdf.cell(0, 30, "Ethnicity Distribution", 0, 1, 'L')
	pdf.set_font('Arial', '', 16)
	pdf.cell(0, 10, 'Ethnic Background of all samples:', 0, 1, 'L')
	
	for key in all_samples_ethnicity: # writes out to PDF of number of samples in each ethnic group
		pdf.set_x(20)
		pdf.cell(0, 8, str(key)+':  '+str(all_samples_ethnicity[key]), 0, 1, 'L')
	
	pdf.cell(0, 10, 'Ethnic Background of removed samples:', 0, 1, 'L')
	for key in ethnicity_removals: # write out to PDF the ethnic backgrounds of the removed samples
		pdf.set_x(20)
		pdf.cell(0, 8, str(key)+':  '+str(ethnicity_removals[key]), 0, 1, 'L')

	# create a files called "samples_to_remove.txt" to be passed in proper format to plink for sample removal
	temp_remove = {x:list(sample_qc_table.loc[sample_qc_table['Sample ID'] == x]['FID']) for x in samples_to_remove}
	for key in temp_remove:
		samples_to_remove_text.write(str(temp_remove[key][0]) + '\t' + str(key) + '\n')
	samples_to_remove_text.flush() # flushes out buffer

	return samples_to_remove_text


def illumina_snp_overview(inputFile, pdf):
	snp_qc_table = pandas.read_table(inputFile, delim_whitespace=True)

	return snps_to_remove
