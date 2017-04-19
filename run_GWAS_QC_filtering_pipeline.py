from chunkypipes.components import *
import sys
import os


class Pipeline(BasePipeline):
	
	def dependencies(self):
		# assuming user as pip installed
		return ['pandas', 'matplotlib', 'fpdf', 'Pillow', 'seaborn', 'pypdf2']

	def description(self):
		return 'Pipeline to perform sample QC (call rate, HWE, Mendelian Error)'

	def configure(self):
			
		return {
			'plink':{
				'path': 'Full path to PLINK executable (must be version >=1.9):'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-sampleTable', required=True, type=str, help="[REQUIRED] Full path to text file of Illumina sample table metrics tab-delimited")
		parser.add_argument('-snpTable', required=True, type=str, help="[REQUIRED] Full path to text file of Illumina SNP table tab-delimited")
		parser.add_argument('-inputPLINK', required=True, type=str, help="Full path to PLINK file to be used in analysis corresponding MAP files or .bim,.fam should be located in same directory (ends in .PED or .BED)")
		parser.add_argument('--arrayType', default='Illumina MEGA', type=str, help='Name of array or chip used for SNPs')
		parser.add_argument('--outDir', default=os.getcwd(), type=str, help='[default:current working directory] Full path to output directory, (note a new directory is made in this directory')
		parser.add_argument('--projectName', default='test', type=str, help="Name of project or owner of project")
		parser.add_argument('--callrate', default=0.97, type=float, help="[default:0.991] minimum call rate to be included in sample set")
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
		parser.add_argument('--chipFailure', default=1, type=int, help='[default:1] Maximum number of sex discrepencies a chip can have before considered failing')

	@staticmethod
	def check_input_format(inputPlinkfile, plink):
		if inputPlinkfile[-4:].lower() == '.ped':
			print "Input .ped, converting to binary"
			convert = subprocess.call(['./'+str(plink), '--file', str(inputPlinkfile[:-4]), '--make-bed', '--out', str(inputPlinkfile[:-4])])

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
		snp_fails_dataframe['F_MISS'] = 'failed '+ str(chrm) + ' SNP call rate threshold: ' + snp_fails_dataframe['F_MISS'].astype(str)
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
		pdf.set_font('Arial', 'B', 14)
		pdf.set_fill_color(200)
		pdf.multi_cell(0, 8, "Total "+str(chrm)+" SNPs analyzed: " +str(total_snps), 1, 'L', True)
		pdf.multi_cell(0, 8, "Total samples considered: " +str(total_samples), 1, 'L', True)
		pdf.multi_cell(0, 8, "Total "+str(chrm)+" SNPs passing call rate threshold:  "+str(total_snps-len(snp_fails)) + '  ' 
				+ '('+str((float(total_snps-len(snp_fails))/float(total_snps))*100)+'%)', 1, 'L', True)
		pdf.set_font('Arial', '', 14)
		pdf.multi_cell(0, 8, "Summary Stats on Original Data:", 1, 1, 'L')
		pdf.set_x(40)
		pdf.multi_cell(0, 8, "Median " + str(chrm) + " missing call rate:  "+ str(stats.median(list(missingness_snp['F_MISS']))*100)+'%', 1, 1, 'L')
		pdf.set_x(40)
		pdf.multi_cell(0, 8, "Mean " + str(chrm) + " missing call rate:  "+ str(stats.mean(list(missingness_snp['F_MISS']))*100)+'%', 1, 1, 'L')
		pdf.set_x(40)	
		pdf.multi_cell(0, 8, "Standard deviation of " + str(chrm) + " missing call rate:  "+ str(stats.stdev(list(missingness_snp['F_MISS']))*100)+'%', 1, 1, 'L')
		pdf.set_x(40)
		pdf.multi_cell(0, 8, "Minimum " + str(chrm) + " missing call rate:  "+ str(min(list(missingness_snp['F_MISS']))*100)+'%', 1, 1, 'L')
		pdf.set_x(40)
		pdf.multi_cell(0, 8, "Maximum " + str(chrm) + " missing call rate:  "+ str(max(list(missingness_snp['F_MISS']))*100)+'%', 1, 1, 'L')
		
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'L')

		return snps_to_remove, remove_reasons

	
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
		pdf_title.cell(0, 10, 'Univerisity of Colorado', 0, 1, 'C')
		pdf_title.cell(0, 10, 'Anschutz Medical Campus', 0, 1, 'C')
		pdf_title.multi_cell(0, 10, '\n\n\n\n', 0, 1, 'C')
		pdf_title.set_font('Arial', 'B', 16)
		pdf_title.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_title.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_title.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_title.image("/home/tonya/Pictures/transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)

		generate_report.explanation_of_deliverables(pdf=pdf_title, params=pipeline_args)

		# write thresholds and parameters to PDF file
		pdf_thresh = FPDF()
		generate_report.thresholds_and_parameters(pdf=pdf_thresh, params=pipeline_args)
		# a list of files to remove once pipeline in finished running, clean up purposes
		stage_for_deletion = []
		
		# write bulk analysis and statistics data to this file
		pdf = FPDF()
		# *****JUST ILLUMINA BASED STATS HERE, NO ACTUAL FILTERING!*****
		# Illumina Threshold Filters, generate stats and create list of samples/snps to remove
		# no actual removal happens here, just list removal and records statistics in PDF
		sample_qc_table, remove_samples_text, stage_for_deletion = generate_report.illumina_sample_overview(inputFile=pipeline_args['sampleTable'], pdf=pdf, callrate=pipeline_args['callrate'], outDir=outdir, cleanup=stage_for_deletion)
		
		snps_to_remove, reasons_snps_fail = generate_illumina_snp_stats.illumina_snp_overview(inputFile=pipeline_args['snpTable'], pdf=pdf, clusterSep=pipeline_args['clusterSep'], aatmean=pipeline_args['AATmean'],
					aatdev=pipeline_args['AATdev'], bbtmean=pipeline_args['BBTmean'], bbtdev=pipeline_args['BBTdev'], aarmean=pipeline_args['AARmean'], abrmean=pipeline_args['ABRmean'],
					bbrmean=pipeline_args['BBRmean'], callrate=pipeline_args['snp_callrate'], outDir=outdir)


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
			Parameter('--remove', outdir+'/males_and_unknowns_by_ped.txt'),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only')
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
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only'),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only')
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
			missingnessSNP=pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.lmiss', missingnessSample=pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.imiss', 
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
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_autosomes.bed')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_autosomes.bim')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_autosomes.fam')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_autosomes.imiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_autosomes.lmiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.bed')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.bim')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.fam')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.lmiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_X_snps_males_and_females.imiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.bed')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.bim')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.fam')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.bed')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.bim')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.fam')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.lmiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_Y_snps_males_only.imiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.lmiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_MT_snps_females_only.imiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.bed')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.bim')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.fam')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.imiss')
		stage_for_deletion.append(pipeline_args['inputPLINK'][:-4]+'_unknown_chr_snps.lmiss')

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
		stage_for_deletion = generate_report.graph_sexcheck(pdf=overall_sex_pdf, sexcheck=pipeline_args['inputPLINK'][:-4]+'_passing_QC.sexcheck', maxF=pipeline_args['maxFemale'], minM=pipeline_args['minMale'], outDir=outdir, cleanup=stage_for_deletion)
		
		pdf_internal_batch = FPDF()
		# checks sex and call rate at the batch level
		stage_for_deletion, failed_chips, batch_summary = generate_report.batch_effects(pdf=pdf_internal_batch, chipFail=pipeline_args['chipFailure'], sexcheck=pipeline_args['inputPLINK'][:-4]+'_passing_QC.sexcheck', missingness=pipeline_args['inputPLINK'][:-4]+'_passing_QC.imiss', outDir=outdir, cleanup=stage_for_deletion)
		



		# --------------------PDF manipulation stuff-----------------------
		pdf_internal_cover_page = FPDF()
		pdf_internal_cover_page.add_page()
		pdf_internal_cover_page.set_margins(20, 10, 20)
		pdf_internal_cover_page.set_font('Arial', 'B', 30)
		pdf_internal_cover_page.cell(0, 30, "Internal QC Report", 0, 1, 'C')
		pdf_internal_cover_page.set_font('Arial', 'B', 20)
		pdf_internal_cover_page.cell(0, 10, 'Univerisity of Colorado', 0, 1, 'C')
		pdf_internal_cover_page.cell(0, 10,'Anschutz Medical Campus', 0, 1, 'C')
		pdf_internal_cover_page.multi_cell(0, 10, '\n\n\n\n', 0, 1, 'C')
		pdf_internal_cover_page.set_font('Arial', 'B', 16)
		pdf_internal_cover_page.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_internal_cover_page.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_internal_cover_page.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_internal_cover_page.image("/home/tonya/Pictures/transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)


		# create cover page for detailed report
		pdf_detailed_cover = FPDF()
		pdf_detailed_cover.add_page()
		pdf_detailed_cover.set_margins(20, 10, 20)
		pdf_detailed_cover.set_font('Arial', 'B', 30)
		pdf_detailed_cover.cell(0, 30, "Detailed QC Report", 0, 1, 'C')
		pdf_detailed_cover.set_font('Arial', 'B', 20)
		pdf_detailed_cover.cell(0, 10, 'Univerisity of Colorado', 0, 1, 'C')
		pdf_detailed_cover.cell(0, 10,'Anschutz Medical Campus', 0, 1, 'C')
		pdf_detailed_cover.multi_cell(0, 10, '\n\n\n\n\n', 0, 1, 'C')
		pdf_detailed_cover.set_font('Arial', 'B', 16)
		pdf_detailed_cover.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_detailed_cover.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_detailed_cover.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_detailed_cover.image("/home/tonya/Pictures/transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)


		# create cover page for glossary report
		pdf_glossary_cover = FPDF()
		pdf_glossary_cover.add_page()
		pdf_glossary_cover.set_margins(20, 10, 20)
		pdf_glossary_cover.set_font('Arial', 'B', 30)
		pdf_glossary_cover.cell(0, 30, "Parameters, Glossary and Pipeline", 0, 1, 'C')
		pdf_glossary_cover.set_font('Arial', 'B', 20)
		pdf_glossary_cover.cell(0, 10, 'Univerisity of Colorado', 0, 1, 'C')
		pdf_glossary_cover.cell(0, 10,'Anschutz Medical Campus', 0, 1, 'C')
		pdf_glossary_cover.multi_cell(0, 10, '\n\n\n\n\n', 0, 1, 'C')
		pdf_glossary_cover.set_font('Arial', 'B', 16)
		pdf_glossary_cover.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_glossary_cover.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_glossary_cover.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')
		pdf_glossary_cover.image("/home/tonya/Pictures/transparent_AMC_cropped.png", x=10, y=245, w=190, h=45)




		pdf_summary_page = FPDF()
		generate_report.overall_main_page_stats(pdf=pdf_summary_page, originalFile=pipeline_args['inputPLINK'][:-4], cleanedFile=pipeline_args['inputPLINK'][:-4]+'_passing_QC')

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
