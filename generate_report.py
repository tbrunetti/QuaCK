from fpdf import FPDF
import re
import pandas
import matplotlib.pyplot as plt
import statistics as stats
import collections
import numpy as np
import seaborn as sns


def thresholds_and_parameters(pdf, params):
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)
	pdf.cell(0, 30, "Parameters and Thresholds", 0, 1, 'L')
	pdf.line(10, 32, 210, 32)
	pdf.set_font('Arial', '', 16)
	for key in params:
		pdf.cell(0, 8, str(key)+':     '+str(params[key]), 0, 1, 'L')


def illumina_sample_overview(inputFile, pdf, callrate):
	print "Running Illumina Sample QC..."
	samples_to_remove_text = open('samples_to_remove.txt', 'w')
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)

	sample_qc_table = pandas.read_table(inputFile)
	total_samples = len(list(sample_qc_table['Sample ID']))
	# retrieve sample ids of those with missing call rage less than call rate parameter provided by user; default = 0.991
	samples_to_remove = list(sample_qc_table[sample_qc_table['Call Rate'] < callrate]['Sample ID'])
	basic_call_stats = [stats.median(sample_qc_table['Call Rate']), stats.mean(sample_qc_table['Call Rate']), stats.stdev(sample_qc_table['Call Rate']), min(sample_qc_table['Call Rate']), max(sample_qc_table['Call Rate'])]
	pdf.cell(0, 30, "Sample Quality Assessment", 0, 1, 'L')
	pdf.line(10, 32, 210, 32)
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


def graph_sexcheck(pdf, sexcheck):
	print "checking sex concordance"
	pdf.add_page()
	pdf.cell(0, 30, "Sex Concordance Check", 0, 1, 'L')
	sex_check_dataframe = pandas.read_table(sexcheck, delim_whitespace=True)
	sorted_sex_check_dataframe = sex_check_dataframe.sort(['F'], ascending=True)
	sorted_sex_check_dataframe['rank'] = list(range(1, len(list(sorted_sex_check_dataframe['FID']))+1))
	
	sample_sex = sns.lmplot(x='rank', y='F', hue='PEDSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 25})
	sns.plt.suptitle('Sex and F coefficient based on pedigree sex data')
	sample_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig('sample_sex.png', bbox_inches='tight')

	imputed_sex = sns.lmplot(x='rank', y='F', hue='SNPSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 25})
	sns.plt.suptitle('Sex and F coefficient based on imputed sex data')
	imputed_sex.set(xlabel='ranked samples', ylabel='F inbreeing coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig('imputed_sex.png', bbox_inches='tight')

	discrepencies_bw_imputed_and_collected = sns.lmplot(x='rank', y='F', hue='STATUS', data=sorted_sex_check_dataframe, fit_reg=False, palette={'OK':'black', 'PROBLEM':'red'}, scatter_kws={"s": 25})
	sns.plt.suptitle('Sex and F coefficient discrpencies based on imputed and pedigree data')
	sns.plt.subplots_adjust(top=.9)
	discrepencies_bw_imputed_and_collected.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig('discrepencies_sex.png', bbox_inches='tight')

	print sorted_sex_check_dataframe


