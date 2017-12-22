from fpdf import FPDF
import re
import pandas
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import statistics as stats
import collections
import numpy as np
import seaborn as sns
import warnings



# this method will be called last so everything can be calculated
# then rearrange pages in PDF so this become page 2
def overall_main_page_stats(pdf, originalFile, cleanedFile, concordance, dupCon, sexCheck):
	num_samples_analyzed = sum(1 for line in open(originalFile+'.fam'))
	num_samples_qc_pass = sum(1 for line in open(cleanedFile+'.fam'))
	num_snps_analyzed = sum(1 for line in open(originalFile+'.bim'))
	num_snps_qc_pass = sum(1 for line in open(cleanedFile+'.bim'))

	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 30)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "QC Summary Statistics", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 12)
	pdf.multi_cell(0, 5, "Listed below are the basic overall summary statistics for this project \
		broken down into three categories: Sample Summary, SNP Summary, and Data Summary. For more detailed \
		information please refer to the subsequent PDFs and pages.  They will provide you with a more granular view of the \
		quality control pipeline that was performed as well as parameters and thresholds that were used.  At the end \
		of the PDF there is also a definitions page that lists what each metric means and how it was calculated.  If you \
		have any questions or concerns please contact the TICR department at the University of Colorado Anschutz Medical \
		Campus at tonya.brunetti@ucdenver.edu. "+'\n\n\n\n\n', 0, 1, 'J')
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	pdf.multi_cell(0, 8, "Sample Summary", 1, 'L', True)
	pdf.set_x(30)
	pdf.set_font('Arial', '', 16)
	pdf.multi_cell(0, 8, "Total samples analyzed:  "+str(num_samples_analyzed), 1, '1', 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total samples passing QC:  "+str(num_samples_qc_pass), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total samples with gender discrepancies:  "+str(sexCheck), 1, 1, 'L')
	pdf.multi_cell(0, 8, "\n\n", 0, 1, 'L')
	pdf.set_font('Arial', 'B', 16)
	pdf.multi_cell(0, 8, "SNP Summary", 1, 'L', True)
	pdf.set_font('Arial', '', 16)
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total SNPs analyzed:  "+str(num_snps_analyzed), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total SNPs passing QC:  "+str(num_snps_qc_pass), 1, 1, 'L')
	pdf.multi_cell(0, 8, "\n\n", 0, 1, 'L')
	pdf.set_font('Arial', 'B', 16)
	pdf.multi_cell(0, 8, "Data Summary", 1, 'L', True)
	pdf.set_font('Arial', '', 16)
	pdf.set_x(30)
	# this is the total number of snps released samples x snps
	pdf.multi_cell(0, 8, "Total genotypes passing QC: "+str(int(num_samples_qc_pass) * int(num_snps_qc_pass)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Average HapMap trio concordance with 1000 Genomes: " +str(concordance)+'%', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Percent duplicate concordance: "+str(dupCon['percent_concordance'])+'%', 1, 1, 'L')



def explanation_of_deliverables(pdf, params):
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Deliverables", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', 'BI', 16)
	pdf.set_x(20)
	pdf.multi_cell(0, 10, 'Raw Plink Files', 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.set_x(25)
	pdf.multi_cell(0, 5, 'Three PLINK files are provided: .BED, .BIM and .FAM files.  These files contain the RAW PRE-FILTERED \
		samples and SNPs from your project.  These files can be used directly with PLINK software.  For more information regarding \
		the PLINK format please refer to the following website https://www.cog-genomics.org/plink/1.9/'+'\n\n', 0, 1, 'J')

	
	pdf.set_font('Arial', 'BI', 16)
	pdf.set_x(20)
	pdf.multi_cell(0, 10, 'Cleaned Plink and VCF Files', 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.set_x(25)
	pdf.multi_cell(0, 5, 'Three PLINK files are provided: .BED, .BIM, and .FAM files.  These files contain the  \
		samples and SNPs that have passed our QC pipeline from your project.  The only samples that are removed from this data set are \
		those that fail missingness thresholds.  All samples failing gender mismatches and ambiguities remain in the cleaned data set despite failing our QC thresholds \
		in case the investigator can resolve these issues from their manifest file.  On the contrary, all SNPs that fail our QC thresholds are removed from the cleaned \
		data set.  All these files can be used directly with PLINK software. \
		For more information regarding the PLINK format please refer to the following website \
		https://www.cog-genomics.org/plink/1.9/  Additionally, we provide the same cleaned data set as a VCF file  ' + '\n\n', 0, 1, 'J')


	pdf.set_font('Arial', 'BI', 16)
	pdf.set_x(20)
	pdf.multi_cell(0, 10, 'PDF Reports', 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.set_x(25)
	pdf.multi_cell(0, 5, 'Three PDF files are provided:', 0, 1, 'J')
	pdf.set_x(35)
	pdf.multi_cell(0, 5, '\n'+params['projectName']+'_final_summary_report.pdf' + '\n' + params['projectName']+'_final_detailed_report.pdf' +'\n' + 
		params['projectName']+'_final_glossary_report.pdf'+'\n\n', 0, 1, 'L') 
	pdf.set_x(25)
	pdf.multi_cell(0, 5, '     The '+params['projectName']+'_final_summary_report.pdf is this current PDF which contains information \
		regarding the explanation of the deliverables and an overall final QC summary page.  This is meant to provide the investigator \
		with a quick overview of the final number of samples, snps, and genotypes that were analyzed and that passed our QC pipeline. '+'\n' +
		'     The ' + params['projectName']+'_final_detailed_report.pdf contains all of the fine details of the samples and snps that pass \
		the recommended Illumina sample and SNP QC, as well as SNPs that pass call rate thresholds and samples sex concordance.  The \
		information and statistics reported are collected sequentially as the data is pushed through the QC pipeline.' +'\n' + '     Finally, \
		the '+ params['projectName']+'_final_glossary_report.pdf provides the investigator with an explanation of the parameters and \
		thresholds that were used in the pipeline.  We have also provided two detailed images of our QC pipeline work flow.  All the details \
		of how parameters and thresholds are calculated are provided in this document. Additionally, the names of the pipeline parameter names \
		are provided in the event that the investigator would like to download and run our QC pipeline on their own copmuter or server under \
		differet parameters and thresholds.'+'\n\n', 0, 1, 'J')


	pdf.set_font('Arial', 'BI', 16)
	pdf.set_x(20)
	pdf.multi_cell(0, 10, 'Text Files', 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.set_x(25)
	pdf.multi_cell(0, 5, 'Nine text files are provided:', 0, 1, 'J')
	pdf.set_x(35)
	pdf.multi_cell(0, 5, '\n'+'snps_failing_QC_details.txt' + '\n' +  'samples_failing_QC_details.txt' +'\n' + 'GenomeStudio_samples_table.txt' +
		'\n' + 'GenomeStudio_SNPs_table.txt' + '\n' + 'GenomeStudio_final_report.txt' + '\n' +  'trio_reports.txt' + 
		'\n' + 'final_report_statistics_per_sample.txt'+ '\n' + 'md5_check_sum.txt' + '\n\n', 0, 1, 'L') 
	pdf.set_x(25)
	pdf.multi_cell(0, 5, '     The snps_failing_QC_details.txt is a tab-delimited file that contains all the SNPs that were removed due to failing at least \
		one QC check.  Each line in the first column represents a single SNP name.  Any susequent columns in the line are reasons why the SNP failed. \
		There may be more than one reason and each reason is a new column.  You will notice a number followed by the reasoning; this number represents \
		the value that the particular SNP was calculated for that parameter.' +'\n' + '     The samples_failing_QC_details.txt is formatted similarly to the snps_failing_QC_details.txt \
		files, except there is an additional column, which is the first column of the files that contains the family ID followed by the second column which \
		contains the sample ID.' +'\n' + '     Alhtough these samples fail QC, only those failing the missingness threshold are actually removed from the analysis.  \
		The GenomeStudio text files are files that contain some infomation that we use in the initial step of our QC pipeline \
		regarding your samples and SNPs.  For more information on what the columns mean please refer to the glossary report PDF.' +'\n' '     We also provide a \
		separate file called trio_reports.txt that gives information on how well the trios performed against 1000 genomes as well as the number of Mendel errors \
		in each family and the concordance rate with 1000 genomes. For trios in which an individual has failed QC, we attempt to run a duo anlaysis assuming a parent \
		and child still remain.  If this is not the case, we list it in the trio_reports.txt as ERROR: NOT FULL DUO OR TRIO. The final_report_statistics_per_sample.txt \
		provides general statistics of the log R ratio and B allele frequencies at the per sample level of samples that pass QC.  Finally, we provide a md5_check_sums.txt file \
		that lists the md5 check sums of all the generated files to ensure files are not altered or corrupted during the transfer process.'+'\n\n', 0, 1, 'J')

	# add disclaimer about storage
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Data Storage and Analysis Disclaimer", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', 'BI', 16)
	pdf.set_x(20)
	pdf.multi_cell(0, 10, 'Data Storage', 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.set_x(25)
	pdf.multi_cell(0, 5, 'Due to the nature of the size of these data sets, we can only keep them on our FTP servers for a limited time before we \
		archive them to save on space.  With that in mind, we would like you to please take a look at the deliverables and files immediately in \
		the event you need access to other files or need something to be re-run.  You will have 3 WEEKS from the date you received the deliverables \
		to download the data and to email us with any questions or concerns or to download data that was not provided to you in the \
		zipped files.  After this 3 week window we will archive all the data and place it on a storage server for 9 months before deleting it from our system.  \
		Please be aware, if you request access to the archived data you will be charged for download time.'+'\n\n', 0, 1, 'L')
	pdf.set_font('Arial', 'BI', 16)
	pdf.set_x(20)
	pdf.multi_cell(0, 10, 'QC Analysis', 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.set_x(25)
	pdf.multi_cell(0, 5, 'The default parameters set on this pipeline have been tested and optimized for studies that have a minimum of 500 samples. \
		If your data has fewer samples than this we cannot guarantee the parameters are the most optimal.  It is possible that some of the threshold set \
		may need to be less or more stringent.  It is up to you to notify us within 3 WEEKS of receiving the reports if you would like to re-run the pipeline under \
		different threholds.  The pipeline is also publically available and is hosted on the following URL: https://github.com/tbrunetti/GWAS_QC_pipeline if you choose to run the pipeline yourself.  Please \
		be aware if you ask us to re-run the pipeline on our server after the 3 week mark, you will be charged for additional time.'+'\n\n\n\n\n', 0, 1, 'L')
	pdf.multi_cell(0, 5,  'Thank you for your understanding and please feel free to contact us with any questions for concerns.', 0, 1, 'L')


def thresholds_and_parameters(pdf, params):
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Parameters and Thresholds", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 16)
	for key in params:
		if key not in ['sampleTable', 'snpTable', 'projectName', 'config', 'inputPLINK', 'outDir', 'arrayType', 'finalReport']:
			pdf.multi_cell(0, 8, str(key)+':     '+str(params[key]), 0, 1, 'L')


def illumina_sample_overview(inputFile, fam, pdf, callrate, outDir, cleanup):
	warnings.simplefilter(action = "ignore", category = FutureWarning)
	print "Running Illumina Sample QC..."
	
	samples_to_remove_text = open(outDir+'/'+'samples_to_remove.txt', 'w')
	pdf.add_page()
	pdf.set_margins(20, 10, 20)

	pdf.set_font('Arial', 'B', 24)

	sample_qc_table = pandas.read_table(inputFile)
	total_samples = len(list(sample_qc_table['Sample ID']))
	# retrieve sample ids of those with missing call rage less than call rate parameter provided by user; default = 0.991
	samples_to_remove = list(sample_qc_table[sample_qc_table['Call Rate'] < callrate]['Sample ID'])
	basic_call_stats = [stats.median(sample_qc_table['Call Rate']), stats.mean(sample_qc_table['Call Rate']), stats.stdev(sample_qc_table['Call Rate']), min(sample_qc_table['Call Rate']), max(sample_qc_table['Call Rate'])]
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Illumina Sample Quality Assessment", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_fill_color(200)
	pdf.set_font('Arial', 'B', 16)
	pdf.multi_cell(0, 8, "Total samples analyzed:  "+str(total_samples), 1, 'L', True)
	pdf.set_fill_color(200)
	pdf.set_font('Arial', 'B', 16)
	pdf.multi_cell(0, 8, "Number of samples passing missing call rate threshold:  "+str(total_samples-len(samples_to_remove)), 1, 'L', True)
	pdf.set_font('Arial', '', 16)
	pdf.set_x(30)
	pdf.multi_cell(0, 8, 'Median call rate:  '+str("%.2f" % round(basic_call_stats[0]*100, 2)) + '%', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Mean call rate:  "+ str("%.2f" % round(basic_call_stats[1]*100, 2)) + '%', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Standard deviation call rate:  "+ str("%.2f" % round(basic_call_stats[2]*100, 2)) + '%', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Minimum call rate:  "+ str("%.2f" % round(basic_call_stats[3]*100, 2)) + '%', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Maximum call rate:  "+ str("%.2f" % round(basic_call_stats[4]*100, 2)) + '%', 1, 1, 'L')
	

	# store batch and chip number of sample that fails missingness threshold
	chips_fail_missingness_check = {}

	famfile = pandas.read_table(fam, delim_whitespace=True, names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'])
	reason_samples_fail = open(outDir + '/' + 'samples_failing_QC_details.txt', 'w')
	# create a files called "samples_to_remove.txt" to be passed in proper format to plink for sample removal
	temp_remove = {x:list(famfile.loc[famfile['IID'] == x]['FID']) for x in samples_to_remove}
	for key in temp_remove:
		samples_to_remove_text.write(str(temp_remove[key][0]) + '\t' + str(key) + '\n')
		reason_samples_fail.write(str(temp_remove[key][0]) + '\t' + str(key) + '\t' + str('failed Illumina Sample Missingness: ') + str(list(sample_qc_table[sample_qc_table['Sample ID'] == key]['Call Rate'])[0]) + '\n')
		# get chips that are failing missingness key is batch list is chip number
		if re.search('([A-Z]*[a-z]*[0-9]*)-DNA_([A-Z]{1}[0-9]{2}).*', key):
			batch_id = re.search('([A-Z]*[a-z]*[0-9]*)-DNA_([A-Z]{1}[0-9]{2}).*', key)
			try:
				chipID = [str(batch_id.group(2))[1:]]
				chips_fail_missingness_check.setdefault(batch_id.group(1), []).extend(chipID)
			except TypeError: # in the event the dictionary is not iterable (list)
				chips_fail_missingness_check[batch_id.group(1)].append(chipID)
		
	samples_to_remove_text.flush() # flushes out buffer
	reason_samples_fail.flush()
	

	def check_GC_callrate(sampleInfo, cleanup):
		warnings.simplefilter(action = "ignore", category = FutureWarning)

		sample_quality_graph = sns.jointplot('Call Rate','p10 GC', data=sampleInfo, kind="reg")
		plt.suptitle('Overall Sample Quality')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'sample_gccallrate.png')
		plt.close()
		pdf.image(outDir+'/'+"sample_gccallrate.png", x=20, y=120, w=170)
		cleanup.append(outDir+'/'+"sample_gccallrate.png")  # puts image in line for deletion; happens after final PDF has been generated
	
	
	check_GC_callrate(sampleInfo=sample_qc_table, cleanup=cleanup)

	return sample_qc_table, samples_to_remove_text, reason_samples_fail, chips_fail_missingness_check, cleanup



def graph_sexcheck(pdf, reason_samples_fail, sexcheck, maxF, minF, outDir, cleanup):
	warnings.simplefilter(action = "ignore", category = FutureWarning)

	print "checking sex concordance"
	pdf.add_page()
	pdf.set_font('Arial', 'B', 24)
	pdf.set_margins(20, 10, 20)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Overall Sex Concordance Check", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	sex_check_dataframe = pandas.read_table(sexcheck, delim_whitespace=True)
	sorted_sex_check_dataframe = sex_check_dataframe.sort_values(['F'], ascending=True)
	sorted_sex_check_dataframe['rank'] = list(range(1, len(list(sorted_sex_check_dataframe['FID']))+1))
	
	sample_sex = sns.lmplot(x='rank', y='F', hue='PEDSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
	plt.axhline(y=float(maxF), linestyle="--")
	plt.axhline(y=float(minF), linestyle="--")
	plt.suptitle('Sex and F coefficient based on pedigree sex data')
	sample_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'sample_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"sample_sex.png", x=20, y=85, w=79, h=85)
	cleanup.append(outDir+'/'+"sample_sex.png")  # puts image in line for deletion; happens after final PDF has been generated

	imputed_sex = sns.lmplot(x='rank', y='F', hue='SNPSEX', data=sorted_sex_check_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
	plt.axhline(y=float(maxF), linestyle="--")
	plt.axhline(y=float(minF), linestyle="--")
	plt.suptitle('Sex and F coefficient based on imputed sex data')
	imputed_sex.set(xlabel='ranked samples', ylabel='F inbreeing coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'imputed_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"imputed_sex.png", x=110, y=85, w=79, h=85)
	cleanup.append(outDir+'/'+"imputed_sex.png")  # puts image in line for deletion; happens after final PDF has been generated

	discrepencies_bw_imputed_and_collected = sns.lmplot(x='rank', y='F', hue='STATUS', data=sorted_sex_check_dataframe, fit_reg=False, palette={'OK':'black', 'PROBLEM':'red'}, scatter_kws={"s": 20})
	plt.axhline(y=float(maxF), linestyle="--")
	plt.axhline(y=float(minF), linestyle="--")
	plt.suptitle('Discrepencies between imputed and pedigree data')
	plt.subplots_adjust(top=.9)
	discrepencies_bw_imputed_and_collected.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'discrepencies_sex.png', bbox_inches='tight')
	plt.close()
	pdf.image(outDir+'/'+"discrepencies_sex.png", x=20, y=190, w=79, h=85)
	cleanup.append(outDir+'/'+"discrepencies_sex.png")  # puts image in line for deletion; happens after final PDF has been generated
	

	# determines which discrepenices are probably human error prone versus sample quality error
	problem_calls_only = sorted_sex_check_dataframe.loc[sorted_sex_check_dataframe['STATUS'].isin(['PROBLEM'])]
	concordant_calls = sorted_sex_check_dataframe.loc[sorted_sex_check_dataframe['STATUS'].isin(['OK'])]
	fixed_sex = list(problem_calls_only[(problem_calls_only['F'] >= minF) | (problem_calls_only['F'] <= maxF)]['IID'])
	indeterminate_sex = list(problem_calls_only[(problem_calls_only['F'] < minF) | (problem_calls_only['F'] > maxF)]['IID'])
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	pdf.multi_cell(0, 10, 'Total Number of Concordant Samples:  ' +  str(len(concordant_calls.index)), 1, 'L', True)
	pdf.multi_cell(0, 10, 'Total Number of Discrepencies:  '+str(len(indeterminate_sex)), 1, 'L', True)
	pdf.set_font('Arial', '', 16)
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of outlier samples with clear gender mismatches:  '+str(len(fixed_sex)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of ambiguous samples:  '+str(len(indeterminate_sex)-len(fixed_sex)), 1, 1, 'L')


	ambiguous_samples = list(set(indeterminate_sex) - set(fixed_sex)) # ambiguous sample IDs onlyq

	with open(reason_samples_fail.name, 'a+') as sex_outliers:
		for sample in fixed_sex:
			sex_outliers.write(str(list(sorted_sex_check_dataframe[sorted_sex_check_dataframe['IID'] == sample]['FID'])[0]) + '\t' +str(sample) +'\t' + 
				'gender_mismatch: '+str(list(sorted_sex_check_dataframe[sorted_sex_check_dataframe['IID'] == sample]['F'])[0]) +'\n')
		for ambiguous in ambiguous_samples :
			sex_outliers.write(str(list(sorted_sex_check_dataframe[sorted_sex_check_dataframe['IID'] == ambiguous]['FID'])[0]) + '\t' +str(ambiguous) +'\t' + 
				'ambiguous_gender: '+str(list(sorted_sex_check_dataframe[sorted_sex_check_dataframe['IID'] == ambiguous]['F'])[0]) +'\n')

		sex_outliers.flush()

	return sex_outliers, str(len(indeterminate_sex)), cleanup

def batch_effects(pdf, chipFail, sexcheck, missingness, chip_missingness_fails, maxF, minF, outDir, cleanup):
	warnings.simplefilter(action = "ignore", category = FutureWarning)

	batch_summary = FPDF()
	# sex concordance between batches
	batch_summary.add_page()
	batch_summary.set_font('Arial', 'B', 30)
	batch_summary.set_margins(20, 10, 20)
	batch_summary.set_x(20)
	batch_summary.multi_cell(0, 30, "Batch Statistics", 0, 1, 'L')
	batch_summary.line(20, 32, 190, 32)
	
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
	
	
	batch_summary.set_font('Arial', 'B', 16)
	batch_summary.set_fill_color(200)
	batch_summary.multi_cell(0, 10, 'Total Number of Batches:  ' +  str(len(batch_sex)), 1, 'L', True)
	

	# missingness data format data for seaborn boxplot/strip plot
	all_batch_callrate = []
	for key in batch_missing:
		callrate_per_batch = [[str(key), str(batch_missing[key][i][0]), str(batch_missing[key][i][1])] for i in range(0, len(batch_missing[key]))]
		all_batch_callrate = all_batch_callrate + callrate_per_batch
	# missingness data read data into pandas dataframe
	missing_call_dataframe = pandas.DataFrame(all_batch_callrate, columns=['batch', 'missing call rate (%)', 'wellID'])
	missing_call_dataframe['missing call rate (%)']=missing_call_dataframe['missing call rate (%)'].astype(float)*100
	missing_genotypes = sns.boxplot(x='missing call rate (%)', y='batch', data=missing_call_dataframe, color=".8")
	missing_genotypes = sns.stripplot(x='missing call rate (%)', y='batch', data=missing_call_dataframe, jitter=True)
	plt.suptitle('Overall missing call rate per sample across batches')
	plt.tight_layout(pad=2, w_pad=2, h_pad=2)
	plt.savefig(outDir+'/'+'missing_call_rate_samples.png')
	plt.close()
	batch_summary.image(outDir+'/'+'missing_call_rate_samples.png', x=10, y=140, w=190, h=150)
	cleanup.append(outDir+'/'+'missing_call_rate_samples.png')  # puts image in line for deletion; happens after final PDF has been generated
	

	# get sample missingness statistics across batches
	batch_call_averages = []
	batch_call_averages_paired = {}
	for batch_name in batch_missing:
		temp = missing_call_dataframe.loc[missing_call_dataframe['batch'].isin([batch_name])]
		batch_call_averages.append(temp['missing call rate (%)'].mean())
		batch_call_averages_paired[batch_name] = temp['missing call rate (%)'].mean() 
	
	# record chip statistics
	total_chips = 0
	total_chips_fail = 0 # fail is when chip has 2 or more sex discrepancies
	failing_chip_IDs = {}

	# outputs graphs and statistics per batch based on sex	
	for key in batch_sex:
		pdf.add_page()
		pdf.set_font('Arial', 'B', 18)
		pdf.set_margins(20, 10, 20)
		pdf.multi_cell(0, 30, "Sample Statistics for Batch ID: "+str(key), 0, 1, 'L')
		pdf.line(20, 30, 190, 30) 
		batch_dataframe = pandas.DataFrame(batch_sex[key], columns=['Discrepencies', 'PEDSEX', 'SNPSEX', 'F', 'well'])
		sorted_sex_batch_dataframe = batch_dataframe.sort_values(['F'], ascending=True)
		sorted_sex_batch_dataframe['rank'] = list(range(1, len(list(sorted_sex_batch_dataframe['well']))+1))
		
		# all sex based batch analysis
		sample_sex = sns.lmplot(x='rank', y='F', hue='PEDSEX', data=sorted_sex_batch_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
		plt.axhline(y=float(maxF), linestyle="--")
		plt.axhline(y=float(minF), linestyle="--")
		plt.suptitle('Sex and F coefficient based on pedigree sex data')
		sample_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'sample_sex'+str(key)+'.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+"sample_sex"+str(key)+'.png', x=20, y=85, w=79, h=85)
		cleanup.append(outDir+'/'+"sample_sex"+str(key)+'.png')  # puts image in line for deletion; happens after final PDF has been generated

		imputed_sex = sns.lmplot(x='rank', y='F', hue='SNPSEX', data=sorted_sex_batch_dataframe, fit_reg=False, palette={0:'black', 1:'pink', 2:'blue'}, scatter_kws={"s": 20})
		plt.axhline(y=float(maxF), linestyle="--")
		plt.axhline(y=float(minF), linestyle="--")
		plt.suptitle('Sex and F coefficient based on imputed sex data')
		imputed_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'imputed_sex'+str(key)+'.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+"imputed_sex"+str(key)+'.png', x=110, y=85, w=79, h=85)
		cleanup.append(outDir+'/'+"imputed_sex"+str(key)+'.png')  # puts image in line for deletion; happens after final PDF has been generated


		discrepencies_sex = sns.lmplot(x='rank', y='F', hue='Discrepencies', data=sorted_sex_batch_dataframe, fit_reg=False, palette={"OK":'black', "PROBLEM":'red'}, scatter_kws={"s": 20})
		plt.axhline(y=float(maxF), linestyle="--")
		plt.axhline(y=float(minF), linestyle="--")
		plt.suptitle('Discrepencies between imputed and pedigree data')
		discrepencies_sex.set(xlabel='ranked samples', ylabel='F inbreeding coefficient')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+'discrepencies_sex'+str(key)+'.png', bbox_inches='tight')
		plt.close()
		pdf.image(outDir+'/'+"discrepencies_sex"+str(key)+'.png', x=20, y=190, w=79, h=85)
		cleanup.append(outDir+'/'+"discrepencies_sex"+str(key)+'.png')  # puts image in line for deletion; happens after final PDF has been generated

		contradictions_headers = sorted_sex_batch_dataframe['Discrepencies'].value_counts().index.tolist()
		contradictions = sorted_sex_batch_dataframe['Discrepencies'].value_counts()
		
		if 'PROBLEM' not in contradictions_headers:
			pdf.set_font('Arial', 'B', 14)
			pdf.multi_cell(0, 8, "Total Samples in Batch:   "+str(len(batch_sex[key])), 0, 1,'L')
			pdf.multi_cell(0, 8, "Percent Sex Concordance in Batch:  100.0%", 0, 1, 'L')
			pdf.multi_cell(0, 8, "Total Samples with Sex Discrepencies:   "+'0', 0, 1, 'L')

			# calculate total number of chips in batch
			all_wells = list(sorted_sex_batch_dataframe['well'])
			batch_chip_total = []
			for i in all_wells:
				batch_chip_total.append(i[1:])
			total_chips = total_chips + len(list(set(batch_chip_total)))

		elif 'OK' not in contradictions_headers:
			pdf.set_font('Arial', 'B', 14)
			pdf.multi_cell(0, 8, "Total Samples in Batch:   "+str(len(batch_sex[key])), 0, 1, 'L')
			pdf.multi_cell(0, 8, "Percent Sex Concordance in Batch:  '0.0%", 0, 1, 'L')
			pdf.multi_cell(0, 8, "Total Samples with Sex Discrepencies:   "+ str(contradictions['PROBLEM']), 0, 1, 'L')
			problem_wells = list(sorted_sex_batch_dataframe[sorted_sex_batch_dataframe['Discrepencies'] == 'PROBLEM']['well'])
			pdf.multi_cell(0, 8, "Wells with Sex Discrepencies:  " + ', '.join(problem_wells))
			
			# calculate total number of chips in batch
			all_wells = list(sorted_sex_batch_dataframe['well'])
			batch_chip_total = []
			for i in all_wells:
				batch_chip_total.append(i[1:])
			total_chips = total_chips + len(list(set(batch_chip_total)))

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

			# count and record number of failing chips
			ordered_freq_chips = problem_columns.most_common()
			while len(ordered_freq_chips) !=0 and ordered_freq_chips[0][1] > chipFail:
				total_chips_fail=total_chips_fail + 1
				if key in failing_chip_IDs:
					failing_chip_IDs[key] = failing_chip_IDs[key] + [str(ordered_freq_chips[0][0])]
					ordered_freq_chips.pop(0)
				else:
					failing_chip_IDs[key] = [str(ordered_freq_chips[0][0])]
					ordered_freq_chips.pop(0)
			
			
		else:
			pdf.set_font('Arial', 'B', 14)
			pdf.multi_cell(0, 8, "Total Samples in Batch:   "+str(len(batch_sex[key])), 0, 1, 'L')
			pdf.multi_cell(0, 8, "Percent Sex Concordance in Batch:  " + str("%.2f" % round((float(contradictions['OK'])/float(len(batch_sex[key])))*100, 2))+'%', 0, 1, 'L')
			pdf.multi_cell(0, 8, "Total Samples with Sex Discrepencies:   "+ str(contradictions['PROBLEM']), 0, 1, 'L')
			problem_wells = list(sorted_sex_batch_dataframe[sorted_sex_batch_dataframe['Discrepencies'] == 'PROBLEM']['well'])
			pdf.multi_cell(0, 8, "Wells with Sex Discrepencies:  " + ', '.join(problem_wells))

			# calculate total number of chips in batch
			all_wells = list(sorted_sex_batch_dataframe['well'])
			batch_chip_total = []
			for i in all_wells:
				batch_chip_total.append(i[1:])
			total_chips = total_chips + len(list(set(batch_chip_total)))

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

			# count and record number of failing chips
			ordered_freq_chips = problem_columns.most_common()
			while len(ordered_freq_chips) != 0 and ordered_freq_chips[0][1] > chipFail:
				total_chips_fail=total_chips_fail + 1
				if key in failing_chip_IDs:
					failing_chip_IDs[key] = failing_chip_IDs[key] + [str(ordered_freq_chips[0][0])]
					ordered_freq_chips.pop(0)
				else:
					failing_chip_IDs[key] = [str(ordered_freq_chips[0][0])]
					ordered_freq_chips.pop(0)
			
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
			plt.title('Distribution of problematic chips', fontweight='bold')
			plt.xlabel("chip number", fontweight='bold')
			plt.ylabel("Frequency", fontweight='bold')
			plt.tight_layout(pad=2, w_pad=2, h_pad=2)
			plt.savefig(outDir+'/'+'problem_chips'+str(key)+'.png', bbox_inches='tight')
			plt.close()
			pdf.image(outDir+'/'+"problem_chips"+str(key)+'.png', x=110, y=230, w=79, h=42)
			cleanup.append(outDir+'/'+'problem_rows'+str(key)+'.png')  # puts image in line for deletion; happens after final PDF has been generated
			cleanup.append(outDir+'/'+"problem_chips"+str(key)+'.png')  # puts image in line for deletion; happens after final PDF has been generated

	# get number of chips failing missingness threshold
	chip_fails_from_missigness = 0
	for key, value in chip_missingness_fails.iteritems():
		chip, num_occ = collections.Counter(value).most_common(1)[0]
		if num_occ > chipFail:
			most_freq_fails = collections.Counter(value).most_common()
			while len(most_freq_fails) != 0 and most_freq_fails[0][1] > chipFail:
				chip_fails_from_missigness = chip_fails_from_missigness +1
				most_freq_fails.pop(0)
		else:
			continue

	
	batch_summary.set_font('Arial', 'B', 16)
	batch_summary.set_fill_color(200)
	batch_summary.multi_cell(0, 10, 'Total Number of Chips:  ' +  str(total_chips), 1, 'L', True)
	batch_summary.set_font('Arial', '', 14)
	batch_summary.set_x(40)
	batch_summary.multi_cell(0, 10, 'Total Number Failing due to Sex: ' + str(total_chips_fail), 1, 1, 'L')
	batch_summary.set_font('Arial', '', 14)
	batch_summary.set_x(40)
	batch_summary.multi_cell(0, 10, 'Total Number Failing due to Missingness: ' + str(chip_fails_from_missigness), 1, 1, 'L')

	batch_summary.set_font('Arial', 'B', 16)
	batch_summary.set_fill_color(200)
	batch_summary.multi_cell(0, 10, 'Batch Sample Missingness Statistics:  ', 1, 'L', True)
	batch_summary.set_font('Arial', '', 14)
	batch_summary.set_x(40)
	batch_summary.multi_cell(0, 10, "Mean sample missingness across all batches: "+str("%.2f" % round(stats.mean(batch_call_averages), 2))+'%', 1, 1, 'L') 
	if len(batch_call_averages) > 1: # std deviation can't be calculated if less than 2 batches exist due to no variance
		batch_summary.set_x(40)
		batch_summary.multi_cell(0, 10, "Standard Deviation in sample missingness across all batches: "+str("%.2f" % round(stats.stdev(batch_call_averages), 2)), 1, 1, 'L')
	else:
		batch_summary.set_x(40)
		batch_summary.multi_cell(0, 10, "Standard Deviation in sample missingness across all batches: 0.00", 1, 1, 'L')
	batch_summary.set_x(40)
	batch_summary.multi_cell(0, 10, "Batch with lowest missingness rate: "+str(min(batch_call_averages_paired, key=batch_call_averages_paired.get))+' ('+str("%.4f" % round(min(batch_call_averages), 4))+'%)', 1, 1, 'L')
	batch_summary.set_x(40)
	batch_summary.multi_cell(0, 10, "Batch with highest missingness rate: "+str(max(batch_call_averages_paired, key=batch_call_averages_paired.get))+' ('+str("%.4f" % round(max(batch_call_averages), 4))+'%)', 1, 1, 'L')


	return cleanup, failing_chip_IDs, batch_summary
