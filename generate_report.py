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
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 30)
	pdf.multi_cell(0, 30, "Parameters and Thresholds", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 16)
	for key in params:
		pdf.multi_cell(0, 8, str(key)+':     '+str(params[key]), 0, 1, 'L')


def illumina_sample_overview(inputFile, pdf, callrate, outDir):
	print "Running Illumina Sample QC..."
	samples_to_remove_text = open(outDir+'/'+'samples_to_remove.txt', 'w')
	pdf.add_page()
	pdf.set_margins(20, 10, 20)

	pdf.set_font('Arial', 'B', 30)

	sample_qc_table = pandas.read_table(inputFile)
	total_samples = len(list(sample_qc_table['Sample ID']))
	# retrieve sample ids of those with missing call rage less than call rate parameter provided by user; default = 0.991
	samples_to_remove = list(sample_qc_table[sample_qc_table['Call Rate'] < callrate]['Sample ID'])
	basic_call_stats = [stats.median(sample_qc_table['Call Rate']), stats.mean(sample_qc_table['Call Rate']), stats.stdev(sample_qc_table['Call Rate']), min(sample_qc_table['Call Rate']), max(sample_qc_table['Call Rate'])]
	pdf.multi_cell(0, 30, "Sample Quality Assessment", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 16)
	pdf.multi_cell(0, 8, "Total samples analyzed:  "+str(total_samples), 0, 1, 'L')
	pdf.multi_cell(0, 8, "Number of samples passing missing call rate threshold:  " + str(total_samples-len(samples_to_remove)), 0, 1, 'L')	
	pdf.multi_cell(0, 8, "Median call rate:  "+ str(basic_call_stats[0]), 0, 1, 'L')
	pdf.multi_cell(0, 8, "Mean call rate:  "+ str(basic_call_stats[1]), 0, 1, 'L')
	pdf.multi_cell(0, 8, "Standard deviation call rate:  "+ str(basic_call_stats[2]), 0, 1, 'L')
	pdf.multi_cell(0, 8, "Minimum call rate:  "+ str(basic_call_stats[3]), 0, 1, 'L')
	pdf.multi_cell(0, 8, "Maximum missing call rate:  "+ str(basic_call_stats[4]), 0, 1, 'L')
	

	# create a files called "samples_to_remove.txt" to be passed in proper format to plink for sample removal
	temp_remove = {x:list(sample_qc_table.loc[sample_qc_table['Sample ID'] == x]['FID']) for x in samples_to_remove}
	for key in temp_remove:
		samples_to_remove_text.write(str(temp_remove[key][0]) + '\t' + str(key) + '\n')
	samples_to_remove_text.flush() # flushes out buffer


	def check_GC_callrate(sampleInfo):
		sample_quality_graph = sns.jointplot('Call Rate','p10 GC', data=sampleInfo, kind="reg")
		sns.plt.suptitle('Overall Sample Quality')
		sns.plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		sns.plt.savefig(outDir+'/'+'sample_gccallrate.png')
		plt.close()
		pdf.image(outDir+'/'+"sample_gccallrate.png", x=20, y=120, w=170)

	
	def check_ethnicity(sampleInfo):
		pdf.add_page()
		pdf.set_margins(20, 10, 20)
		# ethnicity and race breakdown
		sampleInfo['RACE'].replace(np.nan, 'NaN', regex=True, inplace=True) # converts pandas nans into strings (by default they are ints)
		sampleInfo['RACE'].value_counts(dropna=False).plot(kind='pie', autopct='%.2f', fontsize=16)
		plt.axis('equal')
		plt.savefig(outDir+'/'+'ethinic_breakdown.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+'ethinic_breakdown.png', x=20, y=130, w=130)
		
		all_samples_ethnicity = dict(sampleInfo['RACE'].value_counts(dropna=True)) # creates the value counts of ethnicity into dictionary for easy PDF writing key=eth; value=total samples


		store_removal_ethnicity = []
		for i in samples_to_remove:
			 store_removal_ethnicity.append(list(sampleInfo.loc[sampleInfo['Sample ID'] == i]['RACE']))
		store_removal_ethnicity = [store_removal_ethnicity[i][0] for i in range(0, len(store_removal_ethnicity))]
		ethnicity_removals = collections.Counter(store_removal_ethnicity)
		
		pdf.set_font('Arial', 'B', 20)
		pdf.multi_cell(0, 30, "Ethnicity Distribution", 0, 1, 'L')
		pdf.line(20, 32, 190, 32)
		pdf.set_font('Arial', '', 16)
		pdf.multi_cell(0, 10, 'Ethnic Background of all samples:', 0, 1, 'L')
		
		for key in all_samples_ethnicity: # writes out to PDF of number of samples in each ethnic group
			pdf.set_x(20)
			pdf.multi_cell(0, 8, str(key)+':  '+str(all_samples_ethnicity[key]), 0, 1, 'L')
		
		pdf.multi_cell(0, 10, 'Ethnic Background of removed samples:', 0, 1, 'L')
		for key in ethnicity_removals: # write out to PDF the ethnic backgrounds of the removed samples
			pdf.set_x(20)
			pdf.multi_cell(0, 8, str(key)+':  '+str(ethnicity_removals[key]), 0, 1, 'L')		


	check_GC_callrate(sampleInfo=sample_qc_table)
	check_ethnicity(sampleInfo=sample_qc_table)

	return samples_to_remove_text


def graph_sexcheck(pdf, sexcheck, outDir):
	print "checking sex concordance"
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)
	pdf.set_margins(20, 10, 20)
	pdf.multi_cell(0, 30, "Sex Concordance Check", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	sex_check_dataframe = pandas.read_table(sexcheck, delim_whitespace=True)
	sorted_sex_check_dataframe = sex_check_dataframe.sort(['F'], ascending=True)
	sorted_sex_check_dataframe['rank'] = list(range(1, len(list(sorted_sex_check_dataframe['FID']))+1))
	
	sample_sex = sns.lmplot(x='rank', y='F', hue='PEDSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 25})
	sns.plt.suptitle('Sex and F coefficient based on pedigree sex data')
	sample_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'sample_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"sample_sex.png", x=20, y=35, w=79, h=85)

	imputed_sex = sns.lmplot(x='rank', y='F', hue='SNPSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 25})
	sns.plt.suptitle('Sex and F coefficient based on imputed sex data')
	imputed_sex.set(xlabel='ranked samples', ylabel='F inbreeing coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'imputed_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"imputed_sex.png", x=110, y=35, w=79, h=85)

	discrepencies_bw_imputed_and_collected = sns.lmplot(x='rank', y='F', hue='STATUS', data=sorted_sex_check_dataframe, fit_reg=False, palette={'OK':'black', 'PROBLEM':'red'}, scatter_kws={"s": 25})
	sns.plt.suptitle('Discrepencies between imputed and pedigree data')
	sns.plt.subplots_adjust(top=.9)
	discrepencies_bw_imputed_and_collected.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'discrepencies_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"discrepencies_sex.png", x=20, y=130, w=79, h=85)

	print sorted_sex_check_dataframe


