from fpdf import FPDF
import re
import pandas
import matplotlib.pyplot as plt
import statistics as stats
import collections
import numpy as np
import seaborn as sns


# this method will be called last so everything can be calculated
# then rearrange pages in PDF so this become page 2
def overall_main_page_stats(pdf):
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 30)
	pdf.multi_cell(0, 30, "QC Summary Statistics", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 12)
	pdf.multi_cell(0, 8, "Listed below are the basic overall summary statistics for this project \
		broken down into three categories: Sample Summary, SNP Summary, and Data Summary. For more detailed \
		information please refer to the subsequent pages.  They will provide you with a more granular view of the \
		quality control pipeline that was performed as well as parameters and thresholds that were used.  At the end \
		of the PDF there is also a definitions page that lists what each metric means and how it was calculated.  If you \
		have any questions or concerns please contact the TICR department at the University of Colorado Anschutz Medical \
		campus at (303) 724-9918. ", 0, 1, 'L')
	pdf.set_font('Arial', '', 16)
	pdf.multi_cell(0, 8, "Sample Summary", 0, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total Samples analyzed:  ", 0, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total Samples Passing QC:  ", 0, 1, 'L')
	pdf.multi_cell(0, 8, "SNP Summary", 0, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total SNPs analyzed:  ", 0, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total SNPs Passing QC: ", 0, 1, 'L')
	pdf.multi_cell(0, 8, "Data Summary", 0, 1, 'L')
	pdf.set_x(30)
	# this is the total number of snps released samples x snps
	pdf.multi_cell(0, 8, "Total Genotypes Passing QC: ", 0, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Percent Missing Data: ", 0, 1, 'L')



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
	pdf.set_x(40)
	pdf.multi_cell(0, 8, "Median call rate:  "+ str(basic_call_stats[0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 8, "Mean call rate:  "+ str(basic_call_stats[1]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 8, "Standard deviation call rate:  "+ str(basic_call_stats[2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 8, "Minimum call rate:  "+ str(basic_call_stats[3]), 0, 1, 'L')
	pdf.set_x(40)
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
		pdf.multi_cell(0, 30, "Race and Ethnicity Distribution", 0, 1, 'L')
		pdf.line(20, 32, 190, 32)
		pdf.set_font('Arial', '', 16)
		pdf.multi_cell(0, 10, 'Ethnic Background of all samples:', 0, 1, 'L')
		
		for key in all_samples_ethnicity: # writes out to PDF of number of samples in each ethnic group
			pdf.set_x(40)
			pdf.multi_cell(0, 8, str(key)+':  '+str(all_samples_ethnicity[key]), 0, 1, 'L')
		
		pdf.multi_cell(0, 10, 'Ethnic Background of removed samples:', 0, 1, 'L')
		for key in ethnicity_removals: # write out to PDF the ethnic backgrounds of the removed samples
			pdf.set_x(40)
			pdf.multi_cell(0, 8, str(key)+':  '+str(ethnicity_removals[key]), 0, 1, 'L')		


	check_GC_callrate(sampleInfo=sample_qc_table)
	check_ethnicity(sampleInfo=sample_qc_table)

	return samples_to_remove_text

def graph_missingness(sample, snp, pdf, outDir):
	pass;




def graph_sexcheck(pdf, sexcheck, outDir):
	print "checking sex concordance"
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)
	pdf.set_margins(20, 10, 20)
	pdf.multi_cell(0, 30, "Overall Sex Concordance Check", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	sex_check_dataframe = pandas.read_table(sexcheck, delim_whitespace=True)
	sorted_sex_check_dataframe = sex_check_dataframe.sort(['F'], ascending=True)
	sorted_sex_check_dataframe['rank'] = list(range(1, len(list(sorted_sex_check_dataframe['FID']))+1))
	
	sample_sex = sns.lmplot(x='rank', y='F', hue='PEDSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
	sns.plt.suptitle('Sex and F coefficient based on pedigree sex data')
	sample_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'sample_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"sample_sex.png", x=20, y=35, w=79, h=85)

	imputed_sex = sns.lmplot(x='rank', y='F', hue='SNPSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
	sns.plt.suptitle('Sex and F coefficient based on imputed sex data')
	imputed_sex.set(xlabel='ranked samples', ylabel='F inbreeing coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'imputed_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"imputed_sex.png", x=110, y=35, w=79, h=85)

	discrepencies_bw_imputed_and_collected = sns.lmplot(x='rank', y='F', hue='STATUS', data=sorted_sex_check_dataframe, fit_reg=False, palette={'OK':'black', 'PROBLEM':'red'}, scatter_kws={"s": 20})
	sns.plt.suptitle('Discrepencies between imputed and pedigree data')
	sns.plt.subplots_adjust(top=.9)
	discrepencies_bw_imputed_and_collected.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'discrepencies_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"discrepencies_sex.png", x=20, y=130, w=79, h=85)

	

def batch_effects(pdf, sexcheck, missingness, outDir):
	# sex concordance between batches
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)
	pdf.set_margins(20, 10, 20)
	pdf.multi_cell(0, 30, "Batch Statistics", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	
	# sex check and sample missingness by batch and by chip
	sex_check_dataframe = pandas.read_table(sexcheck, delim_whitespace=True)
	sample_missingness_dataframe = pandas.read_table(missingness, delim_whitespace=True)
	batch_sex = {}
	batch_missing = {}
	for batch in list(sex_check_dataframe['IID']):
		if re.search('([A-Z]*[a-z]*[0-9]*)-DNA_([A-Z]{1}[0-9]{2}).*', batch):
			batch_id = re.search('([A-Z]*[a-z]*[0-9]*)-DNA_([A-Z]{1}[0-9]{2}).*', batch)
			if batch_id.group(1) in batch_sex:
				batch_sex[batch_id.group(1)] = batch_sex[batch_id.group(1)] + [list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['STATUS']) + list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['PEDSEX']) + list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['SNPSEX']) + list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['F']) + [batch_id.group(2)]]
				batch_missing[batch_id.group(1)] = batch_missing[batch_id.group(1)] + [list(sample_missingness_dataframe[sample_missingness_dataframe['IID'] == batch]['F_MISS']) + [batch_id.group(2)]]
			else:
				batch_sex[batch_id.group(1)] = [list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['STATUS']) + list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['PEDSEX']) + list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['SNPSEX']) + list(sex_check_dataframe[sex_check_dataframe['IID'] == batch]['F']) + [batch_id.group(2)]]
				batch_missing[batch_id.group(1)] = [list(sample_missingness_dataframe[sample_missingness_dataframe['IID'] == batch]['F_MISS']) + [batch_id.group(2)]]
		else:
			print 'BATCH ID [' + str(batch) + '] not formatted properly!'
	
	
	pdf.set_font('Arial', '', 16)
	pdf.multi_cell(0, 10, 'Total Number of Batches:  ' +  str(len(batch_sex)), 0, 1, 'L')
	# missingness data format data for seaborn boxplot/strip plot
	all_batch_callrate = []
	for key in batch_missing:
		callrate_per_batch = [[str(key), str(batch_missing[key][i][0]), str(batch_missing[key][i][1])] for i in range(0, len(batch_missing[key]))]
		all_batch_callrate = all_batch_callrate + callrate_per_batch
	# missingness data read data into pandas dataframe
	missing_call_dataframe = pandas.DataFrame(all_batch_callrate, columns=['batch', 'missing call rate', 'wellID'])
	missing_call_dataframe['missing call rate']=missing_call_dataframe['missing call rate'].astype(float)
	missing_genotypes = sns.boxplot(x='missing call rate', y='batch', data=missing_call_dataframe, color=".8")
	missing_genotypes = sns.stripplot(x='missing call rate', y='batch', data=missing_call_dataframe, jitter=True)
	sns.plt.suptitle('Overall missing call rate per sample across batches')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'missing_call_rate_samples.png')
	plt.close()
	pdf.image(outDir+'/'+'missing_call_rate_samples.png', x=10, y=130, w=190, h=150)

	# outputs graphs and statistics per batch based on sex	
	for key in batch_sex:
		pdf.add_page()
		pdf.set_font('Arial', 'B', 18)
		pdf.set_margins(20, 10, 20)
		pdf.multi_cell(0, 30, "Sample Statistics for Batch ID: "+str(key), 0, 1, 'L')
		pdf.line(20, 30, 190, 30) 
		batch_dataframe = pandas.DataFrame(batch_sex[key], columns=['Discrepencies', 'PEDSEX', 'SNPSEX', 'F', 'well'])
		sorted_sex_batch_dataframe = batch_dataframe.sort(['F'], ascending=True)
		sorted_sex_batch_dataframe['rank'] = list(range(1, len(list(sorted_sex_batch_dataframe['well']))+1))
		
		# all sex based batch analysis
		sample_sex = sns.lmplot(x='rank', y='F', hue='PEDSEX', data=sorted_sex_batch_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
		sns.plt.suptitle('Sex and F coefficient based on pedigree sex data')
		sample_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'sample_sex'+str(key)+'.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+"sample_sex"+str(key)+'.png', x=20, y=85, w=79, h=85)

		imputed_sex = sns.lmplot(x='rank', y='F', hue='SNPSEX', data=sorted_sex_batch_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
		sns.plt.suptitle('Sex and F coefficient based on imputed sex data')
		imputed_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'imputed_sex'+str(key)+'.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+"imputed_sex"+str(key)+'.png', x=110, y=85, w=79, h=85)

		discrepencies_sex = sns.lmplot(x='rank', y='F', hue='Discrepencies', data=sorted_sex_batch_dataframe, fit_reg=False, palette={"OK":'black', "PROBLEM":'red'}, scatter_kws={"s": 20})
		sns.plt.suptitle('Discrepencies between imputed and pedigree data')
		discrepencies_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'discrepencies_sex'+str(key)+'.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+"discrepencies_sex"+str(key)+'.png', x=20, y=190, w=79, h=85)

		contradictions_headers = sorted_sex_batch_dataframe['Discrepencies'].value_counts().index.tolist()
		contradictions = sorted_sex_batch_dataframe['Discrepencies'].value_counts()
		
		if 'PROBLEM' not in contradictions_headers:
			pdf.set_font('Arial', 'B', 14)
			pdf.multi_cell(0, 8, "Total Samples in Batch:   "+str(len(batch_sex[key])), 0, 1,'L')
			pdf.multi_cell(0, 8, "Percent Sex Concordance in Batch:  100.0%", 0, 1, 'L')
			pdf.multi_cell(0, 8, "Total Samples with Sex Discrepencies:   "+'0', 0, 1, 'L')

		elif 'OK' not in contradictions_headers:
			pdf.set_font('Arial', 'B', 14)
			pdf.multi_cell(0, 8, "Total Samples in Batch:   "+str(len(batch_sex[key])), 0, 1, 'L')
			pdf.multi_cell(0, 8, "Percent Sex Concordance in Batch:  '0.0%", 0, 1, 'L')
			pdf.multi_cell(0, 8, "Total Samples with Sex Discrepencies:   "+ str(contradictions['PROBLEM']), 0, 1, 'L')
			problem_wells = list(sorted_sex_batch_dataframe[sorted_sex_batch_dataframe['Discrepencies'] == 'PROBLEM']['well'])
			pdf.multi_cell(0, 8, "Wells with Sex Discrepencies:  " + ', '.join(problem_wells))
		else:
			pdf.set_font('Arial', 'B', 14)
			pdf.multi_cell(0, 8, "Total Samples in Batch:   "+str(len(batch_sex[key])), 0, 1, 'L')
			pdf.multi_cell(0, 8, "Percent Sex Concordance in Batch:  " + str((float(contradictions['OK'])/float(len(batch_sex[key])))*100)+'%', 0, 1, 'L')
			pdf.multi_cell(0, 8, "Total Samples with Sex Discrepencies:   "+ str(contradictions['PROBLEM']), 0, 1, 'L')
			problem_wells = list(sorted_sex_batch_dataframe[sorted_sex_batch_dataframe['Discrepencies'] == 'PROBLEM']['well'])
			pdf.multi_cell(0, 8, "Wells with Sex Discrepencies:  " + ', '.join(problem_wells))

			rows = []
			columns = []
			for i in problem_wells:
				rows.append(i[0])
				columns.append(i[1:])
			width = 1
			all_rows =['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
			problem_rows = collections.Counter(rows)
			for row in all_rows:
				if row not in problem_rows:
					problem_rows[row] = 0
			
			all_columns=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
			problem_columns = collections.Counter(columns)
			for col in all_columns:
				if col not in problem_columns:
					problem_columns[col] = 0
			
			plt.bar(np.arange(len(problem_rows.keys())), problem_rows.values(), width, edgecolor='black', align='center')
			plt.xticks(np.arange(len(problem_rows.keys())), problem_rows.keys(), fontweight='bold')
			plt.title('Distribution of problematic rows', fontweight='bold')
			plt.xlabel("plate row ID", fontweight='bold')
			plt.ylabel("Frequency", fontweight='bold')
			plt.tight_layout(pad=2, w_pad=2, h_pad=2)
			plt.savefig(outDir+'/'+'problem_rows'+str(key)+'.png', bbox_inches='tight')
			plt.close()
			pdf.image(outDir+'/'+"problem_rows"+str(key)+'.png', x=110, y=190, w=79, h=42)
			plt.bar(np.arange(len(problem_columns.keys())), problem_columns.values(), width, edgecolor='black', align='center')
			plt.xticks(np.arange(len(problem_columns.keys())), problem_columns.keys(), fontweight='bold')
			plt.title('Distribution of problematic columns', fontweight='bold')
			plt.xlabel("plate column number", fontweight='bold')
			plt.ylabel("Frequency", fontweight='bold')
			plt.tight_layout(pad=2, w_pad=2, h_pad=2)
			plt.savefig(outDir+'/'+'problem_columns'+str(key)+'.png', bbox_inches='tight')
			plt.close()
			pdf.image(outDir+'/'+"problem_columns"+str(key)+'.png', x=110, y=230, w=79, h=42)
