from chunkypipes.components import *
import datetime
import sys
import os


class Pipeline(BasePipeline):
	
	def dependencies(self):
		# assuming user as pip installed
		return ['pandas', 'matplotlib', 'fpdf', 'Pillow', 'seaborn', 'statistics', 'pypdf2']

	def description(self):
		return 'Pipeline to perform sample QC (call rate, HWE, Mendelian Error)'

	def configure(self):
			
		return {
			'plink':{
				'path': 'Full path to PLINK executable (must be version >=1.9):'
			},
			'thousand_genomes':{
				'path': 'Full path to 1000 genomes PLINK data ending in .bed for HapMap Trio concordance check:'
			},
			'hapmap_info':{
				'path': 'Full path to tab delimited file containing HapMap FID information (should have 4 headers: IID, SEX, FID, REL):'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-sampleTable', required=True, type=str, help="[REQUIRED] Full path to text file of Illumina sample table metrics tab-delimited")
		parser.add_argument('-snpTable', required=True, type=str, help="[REQUIRED] Full path to text file of Illumina SNP table tab-delimited")
		parser.add_argument('-inputPLINK', required=True, type=str, help="Full path to PLINK file to be used in analysis corresponding MAP files or .bim,.fam should be located in same directory (ends in .PED or .BED)")
		parser.add_argument('--finalReport', default=None, type=str, help='[default:None] Comma separated list. Ex. Full path to tab-delimited genome studio final report with Log R Ratio and B allele frequency column,line header begins on with starting index at 0,the column number where the Sample_ID is located. line0 = line1, line1=line2, etc... \
																			if header columns begin on the 5th line and sample ID is in the 3rd column, then parameter input should be the following: /path/to/file,4,2 ')
		parser.add_argument('--arrayType', default='Illumina MEGA', type=str, help='Name of array or chip used for SNPs')
		parser.add_argument('--outDir', default=os.getcwd(), type=str, help='[default:current working directory] Full path to output directory, do not add last / to path! (note a new directory is made in this directory')
		parser.add_argument('--projectName', default=str(datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")), type=str, help="Name of project or owner of project")
		parser.add_argument('--callrate', default=0.97, type=float, help="[default:0.97] minimum call rate to be included in sample set")
		parser.add_argument('--snp_callrate', default=0.97, type=float, help='[default:0.97] minimum call rate for SNP to be included in autosomal SNP set (anything below this value will be removed')
		parser.add_argument('--clusterSep', default=0.30, type=float, help='[default:0.30] mimimum allowable cluster separation value in order for SNP to be retained (anything equal to or below this value is removed')
		parser.add_argument('--AATmean', default=0.30, help='[default:0.30] maximum allowable AA T mean threshold in order for SNP to be retained (anything above this value is removed')
		parser.add_argument('--AATdev', default=0.06, type=float, help='[default:0.06] maximum allowable AA T dev threshold in order for SNP to be retained (anything above this value is removed)')
		parser.add_argument('--BBTmean', default=0.70, type=float, help='[default:0.70] minimum allowable BB T mean threshold in order for SNP to be retained (anything below this value is removed)')
		parser.add_argument('--BBTdev', default=0.06, type=float, help='[default:0.06] maximum allowable BB T dev threshold in order for SNP to be retained (anything above this value is removed)')
		parser.add_argument('--AARmean', default=0.20, type=float, help='[default:0.20] minimum allowable AA R mean threshold in order for SNP to be retained (anything equal to or below this value is removed)')
		parser.add_argument('--ABRmean', default=0.20, type=float, help='[default:0.20] minimum allowable AB R mean threshold in order for SNP to be retained (anything equal to or below this value is removed)')
		parser.add_argument('--BBRmean', default=0.20, type=float, help='[default:0.20] minimum allowable BB R mean threshold in order for SNP to be retained (anything equal to or below this value is removed)')
		parser.add_argument('--genome_build', default='b37-hg19', type=str, help='[default:b37-hg19], genome build options: b36-hg18, b37-hg19, b38-hg38')
		parser.add_argument('--maxFemale', default=0.20, type=float, help='[default:0.20] F scores below this value will be imputed as female subjects based on X-chromosome imputation')
		parser.add_argument('--minMale', default=0.80, type=float, help='[default:0.80] F scores above this value will be imputed as male subjects based on X-chromosome imputation')
		parser.add_argument('--chipFailure', default=1, type=int, help='[default:1] Maximum number of sex discrepencies or missigness threshold fails a chip can have before considered failing')


	@staticmethod
	def check_input_format(inputPlinkfile, plink):
		if inputPlinkfile[-4:].lower() == '.ped':
			print "Input .ped, converting to binary"
			convert = subprocess.call([str(plink), '--file', str(inputPlinkfile[:-4]), '--make-bed', '--out', str(inputPlinkfile[:-4])])

		elif inputPlinkfile[-4:].lower() == '.bed':
			print "Input seems to follow bed format"
		
		else:
			sys.exit("Error!! Input not recognized, please input .ped or .bed PLINK file" )

	@staticmethod
	def extract_X_boundries(build):
		if build == 'b37-hg19':
			last_bp_head = 2699520
			first_bp_tail = 154931044
			return last_bp_head, first_bp_tail
		elif build == 'b36-hg18':
			last_bp_head = 2709521
			first_bp_tail = 154584237
			return last_bp_head, first_bp_tail
		elif build == 'b38-hg38':
			last_bp_head = 2781479
			first_bp_tail = 155701383
			return last_bp_head, first_bp_tail
		else:
			print "Specified build does not exist, please select different build (see -h for options)"

	@staticmethod
	def separate_sexes(plinkFam, outDir):
		import pandas

		females_and_unknowns = open(outDir+'/females_and_unknowns_by_ped.txt', 'w')
		males_and_unknowns = open(outDir+'/males_and_unknowns_by_ped.txt', 'w')
		unknowns_only = open(outDir+'/unknown_by_ped.txt', 'w')
		sample_info = pandas.read_table(plinkFam, delim_whitespace=True)
		for row in range(0, len(sample_info)):
			print 
			if sample_info.iloc[row, 4]  == 1 or sample_info.iloc[row, 4] == 0: # also place unknowns in this file
				males_and_unknowns.write('\n'+ str(sample_info.iloc[row, 0]) + '\t' +str(sample_info.iloc[row, 1]))
			elif sample_info.iloc[row, 4] == 2 or sample_info.iloc[row, 4] == 0: # also place unknowns in this file
				females_and_unknowns.write('\n'+ str(sample_info.iloc[row, 0]) + '\t' +str(sample_info.iloc[row, 1]))
			elif sample_info.iloc[row, 4] == 0: # unkown sex
				unknowns_only.write('\n'+ str(sample_info.iloc[row, 0]) + '\t' +str(sample_info.iloc[row, 1]))
			else:
				print str(sample_info[row, 1])+":  SKIPPING samples for SNP call rate sex analysis"

		# ensure results are pushed out of buffer
		females_and_unknowns.flush()
		males_and_unknowns.flush()
		unknowns_only.flush()

	@staticmethod
	def call_rate(missingnessSNP, missingnessSample, callrate, snps_to_remove, remove_reasons, pdf, chrm):
		import pandas
		import statistics as stats
		from fpdf import FPDF
		
		missingness_snp = pandas.read_table(missingnessSNP, delim_whitespace=True)
		samples = pandas.read_table(missingnessSample, delim_whitespace=True)
		total_samples = len(list(samples['IID']))
		total_snps = len(list(missingness_snp['SNP']))
		snp_fails_dataframe = missingness_snp[missingness_snp['F_MISS'] > (1-callrate)]
		snp_fails_dataframe['F_MISS'] = 'failed '+ str(chrm) + ' SNP call rate threshold ratio of missing: ' + snp_fails_dataframe['F_MISS'].astype(str)
		snp_fails_dict = snp_fails_dataframe.set_index('SNP')['F_MISS'].to_dict()
		# to make easy to merge with existing dictionary without overwriting values
		snp_fails_dict_appendable = {key:[value] for key, value in snp_fails_dict.iteritems()}
		snp_fails = [key for key in snp_fails_dict] # subtract 1 since parameter is minimum percent called 
		snps_to_remove = snps_to_remove + snp_fails

		# merge dictionary with existing for snp reasons
		for key, value in snp_fails_dict_appendable.iteritems():
			try:
				remove_reasons.setdefault(key, []).extend(value)
			except TypeError:
				remove_reasons[key].append(value)

		# write this information to PDF
		if chrm == 'chr Y':
			pdf.set_margins(20, 10, 20)
			pdf.set_font('Arial', '', 12)
			pdf.set_x(20)
			pdf.multi_cell(0, 5, "This is the table of call rate statistics for " + str(chrm) + " SNPs which passed Illumina recommended sample and SNP QC.  \
				This was calculated by extracting all " + str(chrm) + ' SNPs from ONLY MALE subjects that pass Illumina QC and calculating the call rate \
				exclusively on the ' + str(chrm) + ' subset of SNPs', 0, 1, 'L')
		else:
			pdf.set_margins(20, 10, 20)
			pdf.set_font('Arial', '', 12)
			pdf.set_x(20)
			pdf.multi_cell(0, 5, "This is the table of call rate statistics for " + str(chrm) + " SNPs which passed Illumina recommended sample and SNP QC.  \
				This was calculated by extracting all " + str(chrm) + ' SNPs from all Illumina QC passing samples and calculating the call rate \
				excusively on the ' + str(chrm) + ' subset of SNPs', 0, 1, 'L')
		pdf.set_font('Arial', 'B', 14)
		pdf.set_fill_color(200)
		pdf.multi_cell(0, 8, "Total "+str(chrm)+" SNPs analyzed: " +str(total_snps), 1, 'L', True)
		pdf.multi_cell(0, 8, "Total samples considered: " +str(total_samples), 1, 'L', True)
		pdf.multi_cell(0, 8, "Total "+str(chrm)+" SNPs passing the call rate threshold:  "+str(total_snps-len(snp_fails)) + '  ' 
				+ '('+str("%.2f" % round((float(total_snps-len(snp_fails))/float(total_snps))*100, 2))+'%)', 1, 'L', True)
		pdf.set_font('Arial', '', 14)
		pdf.multi_cell(0, 8, "Summary Stats on Original Data:", 1, 1, 'L')
		pdf.set_x(40)
		pdf.multi_cell(0, 8, "Median " + str(chrm) + " call rate:  "+ str("%.2f" % round(stats.median(list(1-missingness_snp['F_MISS'].astype(float)))*100, 2))+'%', 1, 1, 'L')
		pdf.set_x(40)
		try:
			pdf.multi_cell(0, 8, "Mean " + str(chrm) + " call rate:  "+ str("%.2f" % round(stats.mean(list(1-missingness_snp['F_MISS'].astype(float)))*100, 2))+'%', 1, 1, 'L')
			pdf.set_x(40)
			pdf.multi_cell(0, 8, "Standard deviation of " + str(chrm) + " call rate:  "+ str("%.2f" % round(stats.stdev(list(1-missingness_snp['F_MISS'].astype(float)))*100, 2))+'%', 1, 1, 'L')
			pdf.set_x(40)
			pdf.multi_cell(0, 8, "Minimum " + str(chrm) + " call rate:  "+ str("%.2f" % round(min(list(1-missingness_snp['F_MISS'].astype(float)))*100, 2))+'%', 1, 1, 'L')
			pdf.set_x(40)
			pdf.multi_cell(0, 8, "Maximum " + str(chrm) + " call rate:  "+ str("%.2f" % round(max(list(1-missingness_snp['F_MISS'].astype(float)))*100, 2))+'%', 1, 1, 'L')
		
			pdf.multi_cell(0, 8, '\n', 0, 1, 'L')
		except AttributeError:
			pdf.multi_cell(0, 8, "Mean " + str(chrm) + " call rate:  NaN", 1, 1, 'L')
                        pdf.set_x(40)
                        pdf.multi_cell(0, 8, "Standard deviation of " + str(chrm) + " call rate:  NaN",1, 1, 'L')
                        pdf.set_x(40)
                        pdf.multi_cell(0, 8, "Minimum " + str(chrm) + " call rate:  NaN", 1, 1, 'L')
                        pdf.set_x(40)
                        pdf.multi_cell(0, 8, "Maximum " + str(chrm) + " call rate:  NaN", 1, 1, 'L')

		del missingness_snp
		del samples
		return snps_to_remove, remove_reasons

	

	@staticmethod
	def concordance(qcPassPLINK, plink, hapmap, tgp, outdir):
		import pandas
		import subprocess
		import re
		import numpy as np

		stage_for_deletion = []
		indConc = []

		'''
		TO DO:
		
		make sure qcPassPLINK is not full path or change code to not join paths!!!  Also, should just be plink prefix (no file extension)

		'''

		concordanceResults = open(os.path.join(outdir, 'individal_concordance_reports.txt'), 'w')
		concordanceResults.write('\t'.join(['HapMap_ID', 'overlapping_snps', 'nonmissing_snps', 'concordant_snps', 'concordance_rate']) + '\n')

		samplesToCheck = open(os.path.join(outdir, 'checkConcordance.txt'), 'w')


		hapmaps = pandas.read_table(hapmap, dtype=str) # needs to be tab-delimted! not whitespace due to occurence of "maternal/paternal grand(/mother/father)" in REL
		hapmapsAvail = list(hapmaps['IID'])
			


		fam_file = pandas.read_table(qcPassPLINK + '.fam', delim_whitespace=True, header=None, names=['FID', 'IID', 'PAT', 'MAT', 'GEN', 'AFF'])
		for index,row in fam_file.iterrows():
			if row['IID'].split('_')[-1].upper() in hapmapsAvail:
				samplesToCheck.write('\t'.join([str(row['FID']), str(row['IID'])]) + '\n')
		samplesToCheck.flush()
		samplesToCheck.close()

		with open(samplesToCheck.name, 'r') as individuals:
			# extract each hapmap control as single individual in plink format
			for line in individuals:
				temp = open(os.path.join(outdir, 'temp_conc.txt'), 'w')
				temp.write(str(line))
				temp.flush()
				temp.close()
				subprocess.call([plink, '--bfile', qcPassPLINK, '--keep', temp.name, '--make-bed', '--out', os.path.join(outdir, line.split('\t')[1].rstrip())])

				# updates fam file to match TGP naming convention so merge can be performed
				with open(os.path.join(outdir, line.split('\t')[1].rstrip() + '.fam'), 'r') as update:
					for control in update:
						tempFile = open(os.path.join(outdir, 'tempFile.txt'), 'w')
						control = control.split()
						control[1] = control[1].split('_')[-1].rstrip().upper()
						control[0] = control[1].rstrip().upper()
						tempFile.write('\t'.join(control))
						tempFile.flush()
						tempFile.close()

				# rename tempFile.txt to .fam file of plink since will not modifiy in-place
				os.rename(os.path.join(outdir, 'tempFile.txt'), os.path.join(outdir, line.split('\t')[1].rstrip() + '.fam'))

				# run plink merge with TGP
				subprocess.call([plink, '--bfile', os.path.join(outdir, line.split('\t')[1].rstrip()),
								'--bmerge', tgp, tgp[:-4] + '.bim', tgp[:-4] + '.fam',
								'--merge-mode', '7',
								'--out', os.path.join(outdir, 'indi_concordance_1000genomes')
								])

				# extract pertanent info from generated log file of concordance 
				extract_lines = subprocess.Popen(['tail', '-4', os.path.join(outdir, 'indi_concordance_1000genomes.log')], stdout=subprocess.PIPE)
				get_concordance = subprocess.check_output(['grep', '^[0-9]'], stdin=extract_lines.stdout)
				# regex to sift through the concorance lines from aboveQ
				overlaps = re.search('([0-9]*)\soverlapping\scalls', get_concordance)
				nonmissing = re.search('([0-9]*)\snonmissing', get_concordance)
				concordant = re.search('([0-9]*)\sconcordant', get_concordance)
				concordant_rate = re.search('for\sa\sconcordance\srate\sof\s(0\.[0-9]*)', get_concordance)
				
				# write results to text file for each hapmap control individual	
				concordanceResults.write(str(line.split()[1]) + '\t' + str(overlaps.group(1)) + '\t' + 
					str(nonmissing.group(1)) + '\t' + str(concordant.group(1)) + '\t' + str(concordant_rate.group(1)) + '\n')

				indConc.append(round(float(concordant_rate.group(1)) * 100, 2))

				stage_for_deletion.append(os.path.join(outdir, line.split('\t')[1].rstrip() + '.bed'))
				stage_for_deletion.append(os.path.join(outdir, line.split('\t')[1].rstrip() + '.bim'))
				stage_for_deletion.append(os.path.join(outdir, line.split('\t')[1].rstrip() + '.fam'))
				stage_for_deletion.append(os.path.join(outdir, line.split('\t')[1].rstrip() + '.hh'))
				stage_for_deletion.append(os.path.join(outdir, line.split('\t')[1].rstrip() + '.log'))

		concordanceResults.flush()
		concordanceResults.close()


		return stage_for_deletion, round(np.mean(indConc), 2)


	@staticmethod
	def check_sum(outdir, projectName):
		with open(outdir + '/' + projectName + '/md5_check_sum.txt', 'a+') as md5sum_files:
			for files in os.listdir(outdir + '/' + projectName):
				if files == 'md5_check_sum.txt':
					continue;
				else:
					subprocess.call(['md5sum', outdir+'/'+projectName+'/'+str(files)], stdout=md5sum_files)
	
	
	def run_pipeline(self, pipeline_args, pipeline_config):
		
		import datetime
		import subprocess
		import statistics as stats
		sys.path.append(".")
		import pandas
		import generate_report
		import generate_illumina_snp_stats
		import PyPDF2
		from fpdf import FPDF
		import re
		import numpy as np

		# specifying output location and conflicting project names of files generated	
		try:
			os.stat(pipeline_args['outDir']+'/'+pipeline_args['projectName'])
			sys.exit("project already exists!!")
		except:
			print "Making new directory called "+str(pipeline_args['projectName']) + ' located in ' + str(pipeline_args['outDir'])
			outdir = pipeline_args['outDir']+'/'+pipeline_args['projectName']
			os.mkdir(outdir)


		# create PDF object for output
		pdf_title = FPDF()
		pdf_title.add_page()
		pdf_title.set_margins(20, 10, 20)
		pdf_title.set_font('Arial', 'B', 30)
		pdf_title.cell(0, 30, "QC Summary Report", 0, 1, 'C')
		pdf_title.set_font('Arial', 'B', 20)
		pdf_title.cell(0, 10, 'University of Colorado', 0, 1, 'C')
		pdf_title.cell(0, 10, 'Anschutz Medical Campus', 0, 1, 'C')
		pdf_title.multi_cell(0, 10, '\n\n\n\n', 0, 1, 'C')
		pdf_title.set_font('Arial', 'B', 16)
		pdf_title.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_title.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_title.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_title.image("transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)

		generate_report.explanation_of_deliverables(pdf=pdf_title, params=pipeline_args)

		# write thresholds and parameters to PDF file
		pdf_thresh = FPDF()
		generate_report.thresholds_and_parameters(pdf=pdf_thresh, params=pipeline_args)
		# a list of files to remove once pipeline in finished running, clean up purposes
		stage_for_deletion = []
		
		# write bulk analysis and statistics data to this file
		pdf = FPDF()
		
		# Software initializations
		# initiate PLINK software
		plink_general = Software('plink', pipeline_config['plink']['path'])


		# checks file format of PLINK file, if not in binary converts to binary
		self.check_input_format(
			inputPlinkfile=pipeline_args['inputPLINK'], plink=pipeline_config['plink']['path']
			)

		# makes text file lists of family and sample IDs based on self-identified sex
		self.separate_sexes(
			plinkFam=pipeline_args['inputPLINK'][:-4]+'.fam', outDir=outdir
			)


		# *****JUST ILLUMINA BASED STATS HERE, NO ACTUAL FILTERING!*****
		# Illumina Threshold Filters, generate stats and create list of samples/snps to remove
		# no actual removal happens here, just list removal and records statistics in PDF
		sample_qc_table, remove_samples_text, reason_samples_fail, sample_fail_locations, stage_for_deletion = generate_report.illumina_sample_overview(inputFile=pipeline_args['sampleTable'], fam=pipeline_args['inputPLINK'][:-4]+'.fam', pdf=pdf, callrate=pipeline_args['callrate'], outDir=outdir, cleanup=stage_for_deletion)
		
			
		snps_to_remove, reasons_snps_fail = generate_illumina_snp_stats.illumina_snp_overview(inputFile=pipeline_args['snpTable'], pdf=pdf, clusterSep=pipeline_args['clusterSep'], aatmean=pipeline_args['AATmean'],
					aatdev=pipeline_args['AATdev'], bbtmean=pipeline_args['BBTmean'], bbtdev=pipeline_args['BBTdev'], aarmean=pipeline_args['AARmean'], abrmean=pipeline_args['ABRmean'],
					bbrmean=pipeline_args['BBRmean'], callrate=pipeline_args['snp_callrate'], outDir=outdir)




		# NOTE! If X is already split, it will appear as an error in log file but it will not affect any of the downstream processes
		#plink_general.run(
		#	Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
		#	Parameter('--split-x', self.extract_X_boundries(pipeline_args['genome_build'])[0], self.extract_X_boundries(pipeline_args['genome_build'])[1]),
		#	Parameter('no-fail'),
		#	Parameter('--make-bed')
		#	)
		# remove Illumina sample and SNP initial QC (not including call rate):
		# convert list to temporary file for PLINK
	
		snps_to_remove_illumina = open(outdir+'/snps_to_remove_illumina.txt', 'w')
		snps_to_remove_illumina.write('\n'.join(snps_to_remove))
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
			Parameter('--exclude', os.path.realpath(snps_to_remove_illumina.name)),
			Parameter('--remove', os.path.realpath(remove_samples_text.name)),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC')
			)

		stage_for_deletion.append(outdir+'/snps_to_remove_illumina.txt')
		

		# ----------------------------------------- BEGIN TRIO CHECKS ----------------------------------------------------
		# fam_file and dict_samples are used in both trio and duplication concordance checks
		average_hap_map_concordance = []
		fam_file = pandas.read_table(pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC.fam', delim_whitespace=True, header=None, names=['FID', 'IID', 'PED', 'MAT', 'GEN', 'AFF'])
		dict_samples = dict(zip(fam_file.FID, fam_file.IID))

		# get FID, IID of trio locations:
		get_trios = open(outdir + '/get_trios.txt', 'w')
		update_trio_names = open(outdir + '/update_trio_names.txt', 'w')
		
		trios = [(fid, iid, iid.split('_')[-1]) for fid, iid in dict_samples.iteritems() if iid.split('_')[-1][0:2] == 'NA']
		
		hapmap_info_sheet = pandas.read_table(pipeline_config['hapmap_info']['path'], header=0) # store hapmap info sheet as dataframe

		
		if len(trios) > 0:
			
			trios_run_together_update = {} # key is batch, value is list of lists of FID, IID, NAxxxx, NAxxxx
			trios_run_together = {} # key is batch,value is list of lists of FID and IID
			trio_ids_only = {x[2]:[x[0], x[1]] for x in trios} # will be useful for calculating Mendel errors

			for samples in trios:
				get_trios.write(str(samples[0]) + '\t' + str(samples[1]) + '\n') # PLINK format to extract samples
				update_trio_names.write(str(samples[0]) + '\t' + str(samples[1]) + '\t' + str(samples[2]) + '\t' + str(samples[2]) + '\n') # PLINK format for updating FID/IID
			
				batch_id, remainder = samples[1].split('-')

				if batch_id in trios_run_together:
					trios_run_together_update[batch_id].append([str(samples[0]), str(samples[1]), str(samples[2]), str(samples[2])])
					trios_run_together[batch_id].append([str(samples[0]), str(samples[1])])
				else:
					trios_run_together_update[batch_id] = [[str(samples[0]), str(samples[1]), str(samples[2]), str(samples[2])]]
					trios_run_together[batch_id] = [[str(samples[0]),str(samples[1])]]


			get_trios.flush() # push out buffer contents
			update_trio_names.flush() # push out buffer contents


			# TESTING TRIO #
			trio_rates = open(outdir + '/trio_reports.txt', 'w')
			header = ['sample_1', 'sample_2', 'sample_3', 'total_overlapping_calls', 'total_nonmissing', 'total_concordant', 'percent_concordance', 'total_mendel_errors_in_family']
			trio_rates.write('\t'.join(header) + '\n')


			for key,value in trios_run_together.iteritems():
				temp_trio_file_extract = open(outdir + '/temp_trio_file_extract.txt', 'w')
				temp_trio_file_update = open(outdir + '/temp_trio_file_update.txt', 'w')
				if len(value) == 3 or len(value) ==2:
					values = trios_run_together_update.get(key)
					for trio in value:
						temp_trio_file_extract.write('\t'.join(trio) + '\n') # FID, IID	
					for trio in values:
						temp_trio_file_update.write('\t'.join(trio) + '\n')

					temp_trio_file_extract.flush()
					temp_trio_file_update.flush()


					# makes bed of HapMap trio subsets
					plink_general.run(
						Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
						Parameter('--keep', temp_trio_file_extract.name),
						Parameter('--make-bed'),
						Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_temp')
						)
					
					# updates HapMap trio IDs to be concordant with 1000 genomes subsets
					plink_general.run(
						Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_temp'),
						Parameter('--update-ids', temp_trio_file_update.name),
						Parameter('--make-bed'),
						Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_temp_updated')
						)

					# trio concordance check against 1000 genomes
					# Need to extract just NA followed by number part of trios and rename FID and IID with the NA ID (ONLY TO CHECK FOR 1000 GENOMES CONCORDANCE)
					plink_general.run(
						Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_temp_updated'),
						Parameter('--bmerge', pipeline_config['thousand_genomes']['path'][:-4]),
						Parameter('--merge-mode', '7'),
						Parameter('--out', outdir + '/trio_concordance_1000genomes')
						)

					extract_lines = subprocess.Popen(['tail', '-4', outdir + '/trio_concordance_1000genomes.log'], stdout=subprocess.PIPE)
					get_concordance = subprocess.check_output(['grep', '^[0-9]'], stdin=extract_lines.stdout)
					# regex to sift through the concorance lines from above
					overlaps = re.search('([0-9]*)\soverlapping\scalls', get_concordance)
					nonmissing = re.search('([0-9]*)\snonmissing', get_concordance)
					concordant = re.search('([0-9]*)\sconcordant', get_concordance)
					concordant_rate = re.search('for\sa\sconcordance\srate\sof\s(0\.[0-9]*)', get_concordance)
					

					if len(value) == 3:
						for ids in value:
							trio_rates.write(str(ids[1]) + '\t')
						trio_rates.write(str(overlaps.group(1)) + '\t' + str(nonmissing.group(1)) + '\t' + str(concordant.group(1)) + '\t' + str("{0:.2f}".format(float(concordant_rate.group(1))*100)) + '\t')
					else:
						for ids in value:
							trio_rates.write(str(ids[1]) + '\t')
						trio_rates.write('NA' + '\t')
						trio_rates.write(str(overlaps.group(1)) + '\t' + str(nonmissing.group(1)) + '\t' + str(concordant.group(1)) + '\t' + str("{0:.2f}".format(float(concordant_rate.group(1))*100)) + '\t')

					
					average_hap_map_concordance.append(round(float(concordant_rate.group(1))*100, 2))

					update_child = open(outdir + '/update_child.txt', 'w')
					update_sex = open(outdir + '/update_sex.txt', 'w')
					kinship = {}
					hapmap_genders = {}
					update_info = trios_run_together_update.get(key, 'DNE')
					updated_FID = []
					if update_info != 'DNE':
						for trios in update_info:
							print trios
							try:
								fid = hapmap_info_sheet.loc[hapmap_info_sheet['IID'] == trios[2], 'FID'].iloc[0]
								trios[2] = str(fid)
								print trios
								relationship = hapmap_info_sheet.loc[hapmap_info_sheet['IID'] == trios[3], 'REL'].iloc[0]
								print relationship
								kinship[relationship] = trios[3]
								sex = hapmap_info_sheet.loc[hapmap_info_sheet['IID'] == trios[3], 'SEX'].iloc[0]
								if sex == 'Female':
									hapmap_genders[trios[3]] = 'F'
								elif sex == 'Male':
									hapmap_genders[trios[3]] = 'M'
								updated_FID.append(trios)
							except:
								print('ERROR: IID ' + str(trios[2]) + ' not found!')

						print kinship
						update_child.write(str(fid) + '\t' + str(kinship.get('child')) + '\t' + str(kinship.get('father')) + '\t' + str(kinship.get('mother')) + '\n')
						update_child.flush()
						for indivs in hapmap_genders:
							update_sex.write(str(fid) + '\t' + str(indivs) + '\t' + hapmap_genders.get(indivs) +'\n')
						update_sex.flush()


					mendel_renames = open(outdir + '/mendel_renames.txt', 'w')
					for each_hapmap in updated_FID:
						mendel_renames.write('\t'.join(each_hapmap) + '\n')
					mendel_renames.flush()

					plink_general.run(
						Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_temp'),
						Parameter('--update-ids', mendel_renames.name),
						Parameter('--make-bed'),
						Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_mendel_families')
						)


					plink_general.run(
						Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_mendel_families'),
						Parameter('--update-parents', update_child.name),
						Parameter('--make-bed'),
						Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_mendel_families_child_updated')
						)

					plink_general.run(
						Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_mendel_families_child_updated'),
						Parameter('--update-sex', update_sex.name),
						Parameter('--make-bed'),
						Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated_sex_child')
						)

					
					if len(value) == 3:
						# .lmendel is one line per variant [CHR, SNP, N]
						# .imendel one subsection per nuclear family, each subsection has one line per family member [FID, IID, N]
						# .fmendel is one line per nuclear family [FID, PAT, MAT, CHLD, N]	
						plink_general.run(
							Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated_sex_child'),
							Parameter('--mendel'),
							Parameter('--out', outdir + '/mendel_errors_'+str(key))
							)

						try:
							with open(outdir + '/mendel_errors_'+str(key)+'.fmendel') as errors:
								header = next(errors)
								for line in errors:
									line = line.rstrip().split()
							trio_rates.write(str(line[-1]) + '\n')
						except:
							trio_rates.write('NA trio' + '\n')
					
					else:
						# .lmendel is one line per variant [CHR, SNP, N]
						# .imendel one subsection per nuclear family, each subsection has one line per family member [FID, IID, N]
						# .fmendel is one line per nuclear family [FID, PAT, MAT, CHLD, N]	
						plink_general.run(
							Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated_sex_child'),
							Parameter('--mendel'),
							Parameter('--mendel-duos'),
							Parameter('--out', outdir + '/mendel_errors_'+str(key))
							)
						try:
							with open(outdir + '/mendel_errors_'+str(key)+'.fmendel') as errors:
								header = next(errors)
								for line in errors:
									line = line.rstrip().split()
							trio_rates.write(str(line[-1]) + '\n')
						except:
							trio_rates.write('not a parent-child duo' + '\n')

				else:
					values = trios_run_together_update.get(key)
					trio_rates.write('ERROR: NOT FULL DUO OR TRIO' + '\t' + str(values) + '\n')
					# try running as duo


			del hapmap_info_sheet
			trio_rates.flush()
			
			
			stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated.bed', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated.bim', 
								pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated.fam', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios.bed',
								pipeline_args['inputPLINK'][:-4]+'_hapmap_trios.bim', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios.fam', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios.log',
								pipeline_args['inputPLINK'][:-4]+'_hapmap_trios.nosex', pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated.log',
								pipeline_args['inputPLINK'][:-4]+'_hapmap_trios_updated.nosex', outdir + '/mendel_errors_'+str(key)+'.imendel', outdir + '/mendel_errors_'+str(key)+'.lmendel', 
								outdir + '/mendel_errors_'+str(key)+'.fmendel', outdir + '/mendel_errors_'+str(key)+'.hh', outdir + '/update_child.txt', 
								outdir +'/update_sex.txt', outdir +'/trio_concordance_1000genomes.diff', outdir + '/temp_trio_file_extract.txt',
								outdir + '/temp_trio_file_update.txt'])
			


		# ----------------------------------------- END OF TRIO CHECKS ----------------------------------------------------
		# ------------------------------ BEGIN HAPMAP INDIVIDUAL CONCORDANCE CHECKS ---------------------------------------

		tempFiles, avgIndiConc = self.concordance(qcPassPLINK=pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC', plink=pipeline_config['plink']['path'], 
			hapmap=pipeline_config['hapmap_info']['path'], tgp=pipeline_config['thousand_genomes']['path'], outdir=os.path.join(pipeline_args['outDir'], pipeline_args['projectName']))

		stage_for_deletion.extend(tempFiles)


		# ---------------------------------- END HAPMAP INDIVIDUAL CONCORDANCE CHECKS -------------------------------------------
		# ----------------------------------------- BEGIN DUPLICATION CHECKS ----------------------------------------------------
		all_samples = {}
		for fid, iid in dict_samples.iteritems():
			#all_samples.setdefault(iid.split('_')[-1], []).extend([fid, iid]) 
			get_iid = re.search('WG[0-9]*-DNA_[A-H]{1}[0-9]{1,2}_(.*)', iid)
			all_samples.setdefault(get_iid.group(1), []).extend([fid, iid])
		
		# each item in list contains [fid1, iid1, fid2, iid2] where 1 is the duplicate pair of 2
		duplicate_pairs = [value for key, value in all_samples.iteritems() if len(value) > 2] # only extracts duplicates
		
		if len(duplicate_pairs) > 0: # confirms there are duplicates in the data set
			
			extract_dups_1 = open(outdir + '/duplicates1.txt', 'w') # format is FID IID, this will be file to input into PLINK to extact duplicates
			extract_dups_2 = open(outdir + '/duplicates2.txt', 'w') # format is FID IID, this will be file to input into PLINK to extact duplicates
			change_names = open(outdir + '/duplicate_name_updates.txt', 'w') # change names to match, required in order to get concordance calc
			
			for samples in duplicate_pairs:
				extract_dups_1.write(str(samples[0]) + '\t' + str(samples[1]) + '\n') # input format for PLINK
				extract_dups_2.write(str(samples[2]) + '\t' + str(samples[3]) + '\n') # input format for PLINK
				change_names.write(str(samples[2]) + '\t' + str(samples[3]) + '\t' + str(samples[0]) + '\t' + str(samples[1]) + '\n') # change extract_dups_2 ids, to  match dup1
			
			extract_dups_1.flush()
			extract_dups_2.flush()
			change_names.flush()

			# make plink subset of duplicate 1 samples only
			plink_general.run(
				Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
				Parameter('--keep', extract_dups_1.name),
				Parameter('--make-bed'),
				Parameter('--out', pipeline_args['inputPLINK'][:-4] + '_dup1')
				)
			
			# make plink subset of duplicate 2 samples only
			plink_general.run(
				Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
				Parameter('--keep', extract_dups_2.name),
				Parameter('--make-bed'),
				Parameter('--out', pipeline_args['inputPLINK'][:-4] + '_dup2')
				)

			# change names of duplicate 2 to match duplicate 1
			plink_general.run(
				Parameter('--bfile', pipeline_args['inputPLINK'][:-4] + '_dup2'),
				Parameter('--update-ids', change_names.name),
				Parameter('--make-bed'),
				Parameter('--out', pipeline_args['inputPLINK'][:-4] + '_dup2_names_updated')
				)

			# perform concordance checks
			plink_general.run(
				Parameter('--bfile', pipeline_args['inputPLINK'][:-4] + '_dup2_names_updated'),
				Parameter('--bmerge', pipeline_args['inputPLINK'][:-4] + '_dup1'),
				Parameter('--merge-mode', '7'),
				Parameter('--out', outdir + '/duplicate_concordance')
				)


			extract_lines = subprocess.Popen(['tail', '-4', outdir + '/duplicate_concordance.log'], stdout=subprocess.PIPE)
			get_concordance = subprocess.check_output(['grep', '^[0-9]'], stdin=extract_lines.stdout)
			# regex to sift through the concorance lines from above
			overlaps = re.search('([0-9]*)\soverlapping\scalls', get_concordance)
			nonmissing = re.search('([0-9]*)\snonmissing', get_concordance)
			concordant = re.search('([0-9]*)\sconcordant', get_concordance)
			concordant_rate = re.search('for\sa\sconcordance\srate\sof\s(0\.[0-9]*)', get_concordance)
			# this dictionary is eventually passed to a call in generate_report.py (def overall_main_page_stats)
			duplicate_concordance = {'total_overlapping_calls':overlaps.group(1), 'total_nonmissing': nonmissing.group(1), 'total_concordant':concordant.group(1), 'percent_concordance':str("%.2f" % round(float(concordant_rate.group(1))*100, 2))}


			stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4] + '_dup1.bed', pipeline_args['inputPLINK'][:-4] + '_dup1.bim', pipeline_args['inputPLINK'][:-4] + '_dup1.fam', pipeline_args['inputPLINK'][:-4] + '_dup1.hh'])
			stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4] + '_dup2.bed', pipeline_args['inputPLINK'][:-4] + '_dup2.bim', pipeline_args['inputPLINK'][:-4] + '_dup2.fam', pipeline_args['inputPLINK'][:-4] + '_dup2.hh'])
			stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4] + '_dup2_names_updated.bed', pipeline_args['inputPLINK'][:-4] + '_dup2_names_updated.bim', pipeline_args['inputPLINK'][:-4] + '_dup2_names_updated.fam', pipeline_args['inputPLINK'][:-4] + '_dup2_names_updated.hh'])
			stage_for_deletion.extend([outdir + '/update_trio_names.txt'])
			
		else: # there are no duplicates in the data set
			duplicate_concordance = {'total_overlapping_calls': 0, 'total_nonmissing': 0, 'total_concordant':0, 'percent_concordance':str(0.00)}	


		# ----------------------------------------- END OF DUPLICATION CHECKS ----------------------------------------------------

	

		########## make plink files for missingness RECALCULATED POST-ILLUMINA SNP QC ###########
		
		# all samples autosomal only
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
			Parameter('--chr', '1-22'),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_autosomes')
			)

		# females and males, no unknowns, in BED format extract X snps only
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
			Parameter('--chr', 'X'),
			Parameter('--remove', outdir+'/unknown_by_ped.txt'),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females')
			)
		# make males only plink file in BED format and extract Y snps only
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
			Parameter('--chr', 'Y'),
			Parameter('--remove', outdir+'/females_and_unknowns_by_ped.txt'),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only')
			)
		# make females only plink file in BED FORMAT and extract MT snps only
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
			Parameter('--chr', 'MT'),
			Parameter('--remove', outdir+'/unknown_by_ped.txt'),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed')
			)

		# all samples plink file in BED FORMAT and extract unknown chromosomes (0)
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
			Parameter('--chr', '0'),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps')
			)


		# get missingness statistics exculsively from all samples autosomal chrs only 
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_autosomes'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_autosomes')
			)

		# get missingness statistics exculsively from females and males X chromosome
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females')
			)
		# get missingness statistics exculsively from Y males only
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only')
			)
		# get missingness statistics exclusively from MT females only
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed')
			)

		# get missingness statistics exclusively from uknown chromosomal SNPs
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps')
			)

		callrate_pdf = FPDF()
		callrate_pdf.add_page()
		callrate_pdf.set_margins(20, 10, 20)
		callrate_pdf.set_font('Arial', 'B', 24)
		callrate_pdf.set_x(20)
		callrate_pdf.multi_cell(0, 30, "SNP Call Rates", 0, 1, 'L')
		callrate_pdf.line(20, 32, 190, 32)
		
		# autosomal call rate calculate and publish to callrate_pdf
		snps_to_remove, reasons_snps_fail = self.call_rate(
			missingnessSNP=pipeline_args['inputPLINK'][:-4]+'_autosomes.lmiss', missingnessSample=pipeline_args['inputPLINK'][:-4]+'_autosomes.imiss', 
			callrate=pipeline_args['snp_callrate'], snps_to_remove=snps_to_remove, remove_reasons=reasons_snps_fail,
			pdf=callrate_pdf, chrm='autosomal'
			)

		# chrX call rate calculate and publish to callrate_pdf
		snps_to_remove, reasons_snps_fail = self.call_rate(
			missingnessSNP=pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.lmiss', missingnessSample=pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.imiss', 
			callrate=pipeline_args['snp_callrate'], snps_to_remove=snps_to_remove, remove_reasons=reasons_snps_fail,
			pdf=callrate_pdf, chrm='chr X'
			)

		# chrY call rate calculate and publish to callrate_pdf
		snps_to_remove, reasons_snps_fail = self.call_rate(
			missingnessSNP=pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.lmiss', missingnessSample=pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.imiss', 
			callrate=pipeline_args['snp_callrate'], snps_to_remove=snps_to_remove, remove_reasons=reasons_snps_fail,
			pdf=callrate_pdf, chrm='chr Y'
			)

		# chrMT call rate calculate and publish to callrate_pdf
		snps_to_remove, reasons_snps_fail = self.call_rate(
			missingnessSNP=pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.lmiss', missingnessSample=pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.imiss', 
			callrate=pipeline_args['snp_callrate'], snps_to_remove=snps_to_remove, remove_reasons=reasons_snps_fail,
			pdf=callrate_pdf, chrm='chr MT'
			)

		# chrMT call rate calculate and publish to callrate_pdf
		snps_to_remove, reasons_snps_fail = self.call_rate(
			missingnessSNP=pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.lmiss', missingnessSample=pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.imiss', 
			callrate=pipeline_args['snp_callrate'], snps_to_remove=snps_to_remove, remove_reasons=reasons_snps_fail,
			pdf=callrate_pdf, chrm='chr 0 (no chr ID)'
			)

		# tags these files for removal when the pipeline finishes running
		stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4]+'_autosomes.bed', pipeline_args['inputPLINK'][:-4]+'_autosomes.bim', pipeline_args['inputPLINK'][:-4]+'_autosomes.fam', 
			pipeline_args['inputPLINK'][:-4]+'_autosomes.imiss', pipeline_args['inputPLINK'][:-4]+'_autosomes.lmiss'])
		stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.bed', pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.bim', pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.fam', 
			pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.imiss', pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.lmiss'])
		stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.bed', pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.bim', pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.fam', 
			pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.imiss', pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.lmiss'])
		stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.bed', pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.bim', pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.fam', pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.imiss', 
			pipeline_args['inputPLINK'][:-4]+'_MT_snps_unknowns_removed.lmiss'])
		stage_for_deletion.extend([pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.bed', pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.bim', pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.fam', 
			pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.imiss', pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.lmiss'])


		# remove all SNPs from file not passing Illumina QC and not passing call rate QC
		unique_snps_to_remove = set(snps_to_remove)
		all_snps_removed = open(outdir+'/all_snps_removed.txt', 'w')
		all_snps_removed.write('\n'.join(unique_snps_to_remove))
		all_snps_removed.flush()
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC'),
			Parameter('--exclude', os.path.realpath(all_snps_removed.name)),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_passing_QC')
			)


		snps_failing_QC_details = open(outdir + '/snps_failing_QC_details.txt', 'w')
		for key, value in reasons_snps_fail.iteritems():
			snps_failing_QC_details.write(str(key)+'\t'+'\t'.join(value)+'\n')
		snps_failing_QC_details.flush()

		#-----ROUND 1 SNP QC (Illumina recommended SNP metrics threhold removal)------

		

		
		
		
		#----ROUND 1 SAMPLE QC (sex check and batch effects)------

		# perform sex checks
		
		# check sex
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_QC'),
			Parameter('--check-sex', pipeline_args['maxFemale'], pipeline_args['minMale']),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_passing_QC')
			)

		# plink missingness calculations primarily for sample missingness check
		print "Running PLINK sample and snp missingness"
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_QC'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_passing_QC')
			)

		overall_sex_pdf = FPDF()
		# not imbedded in Illumina Sample check because uses own input files
		reason_samples_fail, total_discrepancies, stage_for_deletion = generate_report.graph_sexcheck(pdf=overall_sex_pdf, reason_samples_fail=reason_samples_fail, sexcheck=pipeline_args['inputPLINK'][:-4]+'_passing_QC.sexcheck', maxF=pipeline_args['maxFemale'], minF=pipeline_args['minMale'], outDir=outdir, cleanup=stage_for_deletion)
		
		pdf_internal_batch = FPDF()
		# checks sex and call rate at the batch level
		stage_for_deletion, failed_chips, batch_summary = generate_report.batch_effects(pdf=pdf_internal_batch, chipFail=pipeline_args['chipFailure'], sexcheck=pipeline_args['inputPLINK'][:-4]+'_passing_QC.sexcheck', missingness=pipeline_args['inputPLINK'][:-4]+'_passing_QC.imiss', 
			chip_missingness_fails=sample_fail_locations, maxF=pipeline_args['maxFemale'], minF=pipeline_args['minMale'], outDir=outdir, cleanup=stage_for_deletion)
		



		# --------------------PDF manipulation stuff-----------------------
		pdf_internal_cover_page = FPDF()
		pdf_internal_cover_page.add_page()
		pdf_internal_cover_page.set_margins(20, 10, 20)
		pdf_internal_cover_page.set_font('Arial', 'B', 30)
		pdf_internal_cover_page.cell(0, 30, "Internal QC Report", 0, 1, 'C')
		pdf_internal_cover_page.set_font('Arial', 'B', 20)
		pdf_internal_cover_page.cell(0, 10, 'University of Colorado', 0, 1, 'C')
		pdf_internal_cover_page.cell(0, 10,'Anschutz Medical Campus', 0, 1, 'C')
		pdf_internal_cover_page.multi_cell(0, 10, '\n\n\n\n', 0, 1, 'C')
		pdf_internal_cover_page.set_font('Arial', 'B', 16)
		pdf_internal_cover_page.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_internal_cover_page.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_internal_cover_page.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_internal_cover_page.image("transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)


		# create cover page for detailed report
		pdf_detailed_cover = FPDF()
		pdf_detailed_cover.add_page()
		pdf_detailed_cover.set_margins(20, 10, 20)
		pdf_detailed_cover.set_font('Arial', 'B', 30)
		pdf_detailed_cover.cell(0, 30, "Detailed QC Report", 0, 1, 'C')
		pdf_detailed_cover.set_font('Arial', 'B', 20)
		pdf_detailed_cover.cell(0, 10, 'University of Colorado', 0, 1, 'C')
		pdf_detailed_cover.cell(0, 10,'Anschutz Medical Campus', 0, 1, 'C')
		pdf_detailed_cover.multi_cell(0, 10, '\n\n\n\n\n', 0, 1, 'C')
		pdf_detailed_cover.set_font('Arial', 'B', 16)
		pdf_detailed_cover.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_detailed_cover.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_detailed_cover.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_detailed_cover.image("transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)


		# create cover page for glossary report
		pdf_glossary_cover = FPDF()
		pdf_glossary_cover.add_page()
		pdf_glossary_cover.set_margins(20, 10, 20)
		pdf_glossary_cover.set_font('Arial', 'B', 30)
		pdf_glossary_cover.cell(0, 30, "Parameters, Glossary and Pipeline", 0, 1, 'C')
		pdf_glossary_cover.set_font('Arial', 'B', 20)
		pdf_glossary_cover.cell(0, 10, 'University of Colorado', 0, 1, 'C')
		pdf_glossary_cover.cell(0, 10,'Anschutz Medical Campus', 0, 1, 'C')
		pdf_glossary_cover.multi_cell(0, 10, '\n\n\n\n\n', 0, 1, 'C')
		pdf_glossary_cover.set_font('Arial', 'B', 16)
		pdf_glossary_cover.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_glossary_cover.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_glossary_cover.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_glossary_cover.image("transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)




		pdf_summary_page = FPDF()
		generate_report.overall_main_page_stats(pdf=pdf_summary_page, originalFile=pipeline_args['inputPLINK'][:-4], cleanedFile=pipeline_args['inputPLINK'][:-4]+'_passing_QC', concordance=avgIndiConc, dupCon=duplicate_concordance, sexCheck=total_discrepancies)

		pdf_title.output(outdir + '/'+pipeline_args['projectName']+'_cover_page.pdf', 'F')
		pdf.output(outdir + '/'+pipeline_args['projectName']+'_bulk_data.pdf', 'F')
		pdf_summary_page.output(outdir + '/'+pipeline_args['projectName']+'_summary_page.pdf', 'F')
		pdf_thresh.output(outdir + '/'+pipeline_args['projectName']+'_thresholds.pdf', 'F')
		pdf_internal_batch.output(outdir +'/'+pipeline_args['projectName']+'_internal_batch.pdf', 'F')
		callrate_pdf.output(outdir +'/'+pipeline_args['projectName']+'_non_auto_callrates.pdf', 'F')
		overall_sex_pdf.output(outdir +'/'+pipeline_args['projectName']+'_overall_sex_concordance.pdf', 'F')
		batch_summary.output(outdir + '/' +pipeline_args['projectName']+'_batch_summary.pdf', 'F')
		pdf_internal_cover_page.output(outdir + '/' + pipeline_args['projectName']+'_internal_cover.pdf', 'F')
		pdf_detailed_cover.output(outdir + '/' + pipeline_args['projectName']+'_detailed_cover.pdf', 'F')
		pdf_glossary_cover.output(outdir + '/' + pipeline_args['projectName']+'_glossary_cover.pdf', 'F')

		# create PDF merge objects and write final PDF as project name with '_final_*.pdf' as suffix
		pdf_merger_summary = PyPDF2.PdfFileMerger()
		pdf_merger_detailed = PyPDF2.PdfFileMerger()
		pdf_merger_glossary = PyPDF2.PdfFileMerger()
		pdf_merger_internal = PyPDF2.PdfFileMerger()
		
		# assign pages to approriate PDF merger object in order of page output
		pdf_merger_summary.append(outdir + '/'+pipeline_args['projectName']+'_cover_page.pdf')
		pdf_merger_summary.append(outdir + '/'+pipeline_args['projectName']+'_summary_page.pdf')
		pdf_merger_detailed.append(outdir + '/' + pipeline_args['projectName']+'_detailed_cover.pdf')
		pdf_merger_detailed.append(outdir + '/'+pipeline_args['projectName']+'_bulk_data.pdf')
		pdf_merger_detailed.append(outdir +'/'+pipeline_args['projectName']+'_non_auto_callrates.pdf')
		pdf_merger_detailed.append(outdir + '/'+pipeline_args['projectName']+'_overall_sex_concordance.pdf')
		pdf_merger_glossary.append(outdir + '/' + pipeline_args['projectName']+'_glossary_cover.pdf')
		pdf_merger_glossary.append(outdir + '/'+pipeline_args['projectName']+'_thresholds.pdf')
		pdf_merger_glossary.append('Parameter_Definitions.pdf')
		pdf_merger_glossary.append('pre-QC-initialization-workflow.pdf')
		pdf_merger_glossary.append('QC-pipeline-workflow.pdf')
		pdf_merger_internal.append(outdir + '/' + pipeline_args['projectName']+'_internal_cover.pdf')
		pdf_merger_internal.append(outdir + '/'+ pipeline_args['projectName']+'_batch_summary.pdf')
		pdf_merger_internal.append(outdir + '/'+pipeline_args['projectName']+'_internal_batch.pdf' )
		

		# write out final reports to approriate PDF files
		pdf_merger_summary.write(outdir + '/'+pipeline_args['projectName']+'_final_summary_report.pdf')
		pdf_merger_detailed.write(outdir + '/'+pipeline_args['projectName']+'_final_detailed_report.pdf')
		pdf_merger_glossary.write(outdir + '/'+pipeline_args['projectName']+'_final_glossary_report.pdf')
		pdf_merger_internal.write(outdir + '/'+pipeline_args['projectName']+'_final_internal_report.pdf')

		# close pdf merger objects
		pdf_merger_summary.close()
		pdf_merger_detailed.close()
		pdf_merger_glossary.close()
		pdf_merger_internal.close()

		
		# remove extraneous intermediate files
		stage_for_deletion.append(outdir + '/samples_to_remove.txt')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_cover_page.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_bulk_data.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_summary_page.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_thresholds.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_internal_batch.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_non_auto_callrates.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_overall_sex_concordance.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_batch_summary.pdf')
		stage_for_deletion.append(outdir + '/' + pipeline_args['projectName']+'_internal_cover.pdf')
		stage_for_deletion.append(outdir + '/' + pipeline_args['projectName']+'_detailed_cover.pdf')
		stage_for_deletion.append(outdir + '/' + pipeline_args['projectName']+'_glossary_cover.pdf')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC.bed')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC.bim')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC.fam')
		stage_for_deletion.append(outdir + '/unknown_by_ped.txt')
		stage_for_deletion.append(outdir + '/females_and_unknowns_by_ped.txt')
		stage_for_deletion.append(outdir + '/males_and_unknowns_by_ped.txt')
		stage_for_deletion.append(outdir + '/all_snps_removed.txt')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'*.hh')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'*.lmiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'*.imiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'*.sexcheck')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_passing_Illumina_sample_SNP_QC.*')

		
		# actually remove files in stage_for_deletion
		print '\n\n' + "Cleaning up project directory"
		for files in stage_for_deletion:
			subprocess.call(['rm', '-rf', files])
		
		
		if pipeline_args['finalReport'] != None:
			final_report_stats = open(outdir + '/final_report_statistics_per_sample.txt', 'w')
			final_report_stats.write('\t'.join(['Sample_ID', 'median_LLR', 'mean_LRR', 'std_LRR', 'median_BAF', 'mean_BAF', 'std_BAF', 'max_LLR', 'min_LLR', 'max_BAF', 'min_BAF']) + '\n')
			samples = {}
			path, header, sampleID = pipeline_args['finalReport'].split(',')
			with open(path, 'r') as report:
				f = open(outdir + '/header_stuff.txt', 'w')
				for _ in xrange(int(header)):
					f.write(next(report))
				true_header = next(report)
				for line in report:
					if outdir + '/' + str(line.rstrip().split('\t')[int(sampleID)]) + '.txt' in samples:
						f.write(line + '\n')
					else:
						f.close()
						f = open(outdir + '/' + str(line.rstrip().split('\t')[int(sampleID)]) + '.txt', 'w')
						f.write(true_header)
						samples[f.name] = 1
				f.close()

			
			failing_snp_names = [key for key in reasons_snps_fail]
			for key in samples:
				sample_subset_raw = pandas.read_table(key)
				sample_subset = sample_subset_raw[~sample_subset_raw['SNP Name'].isin(failing_snp_names)]
				del sample_subset_raw
				try:
					final_report_stats.write(str(key.split('/')[-1][:-4]) + '\t' + str(sample_subset['Log R Ratio'].median(skipna=True)) + '\t' + str(sample_subset['Log R Ratio'].mean(skipna=True)) +'\t' +
						str(sample_subset['Log R Ratio'].std(skipna=True)) + '\t' + str(sample_subset['B Allele Freq'].median(skipna=True)) + '\t' + str(sample_subset['B Allele Freq'].mean(skipna=True)) + '\t' +
						str(sample_subset['B Allele Freq'].std(skipna=True)) + '\t' + str(sample_subset['Log R Ratio'].max(skipna=True)) +'\t' + str(sample_subset['Log R Ratio'].min(skipna=True)) + '\t' +
						str(sample_subset['B Allele Freq'].max(skipna=True)) + '\t' + str(sample_subset['B Allele Freq'].min(skipna=True)) +'\n')
					del sample_subset
				except TypeError:
					print 'Wrong type at key ' + str(key)
				except MemoryError:
					print 'Out of Memory at key ' + str(key)

			final_report_stats.close()	

			for key in samples:
				subprocess.call(['rm', '-rf', key])
		
		# creates a cleaned VCF in addition to the cleaned PLINK file
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_passing_QC'),
			Parameter('--recode', 'vcf-iid'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_passing_QC')
			)
	

		self.check_sum(outdir=pipeline_args['outDir'], projectName=pipeline_args['projectName'])
