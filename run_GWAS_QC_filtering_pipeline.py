from chunkypipes.components import *
import sys
import os
import pandas
import datetime
from fpdf import FPDF
sys.path.append(".")
import generate_report
import generate_illumina_snp_stats
import PyPDF2

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
		parser.add_argument('-sampleTable', required=True, help="[REQUIRED] Full path to text file of Illumina sample table metrics tab-delimited")
		parser.add_argument('-snpTable', required=True, help="[REQUIRED] Full path to text file of Illumina SNP table tab-delimited")
		parser.add_argument('-inputPLINK', required=True, help="Full path to PLINK file to be used in analysis corresponding MAP files or .bim,.fam should be located in same directory (ends in .PED or .BED)")
		parser.add_argument('--arrayType', default='MEGAex with Custom Content', help='Name of array or chip used for SNPs')
		parser.add_argument('--outDir', default=os.getcwd(), help='[default:current working directory] Full path to output directory, (note a new directory is made in this directory')
		parser.add_argument('--projectName', default='test', help="Name of project or owner of project")
		parser.add_argument('--callrate', default=0.97, help="[default:0.991] minimum call rate to be included in sample set")
		parser.add_argument('--snp_callrate', default=0.97, help='[default:0.97] minimum call rate for SNP to be included in autosomal SNP set (anything below this value will be removed')
		parser.add_argument('--clusterSep', default=0.30, help='[default:0.30] mimimum allowable cluster separation value in order for SNP to be retained (anything equal to or below this value is removed')
		parser.add_argument('--AATmean', default=0.30, help='[default:0.30] maximum allowable AA T mean threshold in order for SNP to be retained (anything above this value is removed')
		parser.add_argument('--AATdev', default=0.06, help='[default:0.06] maximum allowable AA T dev threshold in order for SNP to be retained (anything above this value is removed)')
		parser.add_argument('--BBTmean', default=0.70, help='[default:0.70] minimum allowable BB T mean threshold in order for SNP to be retained (anything below this value is removed)')
		parser.add_argument('--BBTdev', default=0.06, help='[default:0.06] maximum allowable BB T dev threshold in order for SNP to be retained (anything above this value is removed)')
		parser.add_argument('--AARmean', default=0.20, help='[default:0.20] minimum allowable AA R mean threshold in order for SNP to be retained (anything equal to or below this value is removed)')
		parser.add_argument('--ABRmean', default=0.20, help='[default:0.20] minimum allowable AB R mean threshold in order for SNP to be retained (anything equal to or below this value is removed)')
		parser.add_argument('--BBRmean', default=0.20, help='[default:0.20] minimum allowable BB R mean threshold in order for SNP to be retained (anything equal to or below this value is removed)')
		parser.add_argument('--genome_build', default='b37-hg19', help='[default:b37-hg19], genome build options: b36-hg18, b37-hg19, b38-hg38')
		parser.add_argument('--maxFemale', default=0.20, help='[default:0.20] F scores below this value will be imputed as female subjects based on X-chromosome imputation')
		parser.add_argument('--minMale', default=0.80, help='[default:0.80] F scores above this value will be imputed as male subjects based on X-chromosome imputation')
	

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

	
	def run_pipeline(self, pipeline_args, pipeline_config):
	
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
		pdf_title.cell(0, 30, "QC Report Generated by TICR", 0, 1, 'C')
		pdf_title.set_font('Arial', 'B', 20)
		pdf_title.cell(0, 35, 'Anschutz Medical Campus', 0, 1, 'C')
		pdf_title.cell(0, 10, 'Univerisity of Colorado', 0, 1, 'C')
		pdf_title.set_font('Arial', 'B', 16)
		pdf_title.cell(0, 10, 'Project:  '+ str(pipeline_args['projectName']), 0, 1, 'C')
		pdf_title.cell(0, 10, 'Array/Chip:  '+ str(pipeline_args['arrayType']), 0, 1, 'C')
		pdf_title.cell(0, 10, 'Date:  '+str(datetime.date.today()), 0, 1, 'C')

		
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
		remove_samples, stage_for_deletion = generate_report.illumina_sample_overview(inputFile=pipeline_args['sampleTable'], pdf=pdf, callrate=pipeline_args['callrate'], outDir=outdir, cleanup=stage_for_deletion)
		
		illumina_snps_to_remove = generate_illumina_snp_stats.illumina_snp_overview(inputFile=pipeline_args['snpTable'], pdf=pdf, clusterSep=pipeline_args['clusterSep'], aatmean=pipeline_args['AATmean'],
					aatdev=pipeline_args['AATdev'], bbtmean=pipeline_args['BBTmean'], bbtdev=pipeline_args['BBTdev'], aarmean=pipeline_args['AARmean'], abrmean=pipeline_args['ABRmean'],
					bbrmean=pipeline_args['BBRmean'], callrate=pipeline_args['snp_callrate'], outDir=outdir)


		# TO DO:
		# use plink --remove to remove the samples that fail QC in remove_sample file
		
		# initiate PLINK software
		plink_general = Software('plink', pipeline_config['plink']['path'])

		# checks file format of PLINK file, if not in binary converts to binary
		self.check_input_format(
			inputPlinkfile=pipeline_args['inputPLINK'], plink=pipeline_config['plink']['path']
			)

		#-----ROUND 1 SNP QC (Illumina recommended SNP metrics threhold removal)------

		# remove Illumina sample and SNP initial QC:
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
			Parameter('--exclude', os.path.realpath(illumina_snps_to_remove.name)),
			Parameter('--remove', os.path.realpath(remove_samples.name)),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4]+'_passing_QC')
			)
		
		
		#----ROUND 1 SAMPLE QC (sex check and batch effects)------

		# perform sex checks
		# NOTE! If X is already split, it will appear as an error in log file but it will not affect any of the downstream processes
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
			Parameter('--split-x', self.extract_X_boundries(pipeline_args['genome_build'])[0], self.extract_X_boundries(pipeline_args['genome_build'])[1]),
			Parameter('no-fail'),
			Parameter('--make-bed'),
			)
		# check sex
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
			Parameter('--check-sex', pipeline_args['maxFemale'], pipeline_args['minMale']),
			Parameter('--out', pipeline_args['inputPLINK'][:-4])
			)

		# plink missingness calculations primarily for sample missingness check
		print "Running PLINK sample and snp missingness"
		plink_general.run(
			Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
			Parameter('--missing'),
			Parameter('--out', pipeline_args['inputPLINK'][:-4])
			)

		# not imbedded in Illumina Sample check because uses own input files
		stage_for_deletion = generate_report.graph_sexcheck(pdf=pdf, sexcheck=pipeline_args['inputPLINK'][:-4]+'.sexcheck', outDir=outdir, cleanup=stage_for_deletion)
		# checks sex and call rate at the batch level
		stage_for_deletion = generate_report.batch_effects(pdf=pdf, sexcheck=pipeline_args['inputPLINK'][:-4]+'.sexcheck', missingness=pipeline_args['inputPLINK'][:-4]+'.imiss', outDir=outdir, cleanup=stage_for_deletion)
		

	
		#----ROUND 2 SNP QC-------

		# plink missingness calculations
		#print "Running PLINK sample and snp missingness"
		#plink_general.run(
		#	Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
		#	Parameter('--missing'),
		#	Parameter('--out', pipeline_args['inputPLINK'][:-4])
		#	)


		#----ROUND 2 SAMPLE QC------
	



		# PDF manipulation stuff
		pdf_summary_page = FPDF()
		generate_report.overall_main_page_stats(pdf=pdf_summary_page, originalFile=pipeline_args['inputPLINK'][:-4], cleanedFile=pipeline_args['inputPLINK'][:-4]+'_passing_QC')

		pdf_title.output(outdir + '/'+pipeline_args['projectName']+'_cover_page.pdf', 'F')
		pdf.output(outdir + '/'+pipeline_args['projectName']+'_bulk_data.pdf', 'F')
		pdf_summary_page.output(outdir + '/'+pipeline_args['projectName']+'_summary_page.pdf', 'F')
		pdf_thresh.output(outdir + '/'+pipeline_args['projectName']+'_thresholds.pdf', 'F')
		# create PDF merge object and write final PDF as project name with '_final_report.pdf' as suffix
		pdf_merger = PyPDF2.PdfFileMerger()
		pdf_merger.append(outdir + '/'+pipeline_args['projectName']+'_cover_page.pdf')
		pdf_merger.append(outdir + '/'+pipeline_args['projectName']+'_summary_page.pdf')
		pdf_merger.append(outdir + '/'+pipeline_args['projectName']+'_bulk_data.pdf')
		pdf_merger.append(outdir + '/'+pipeline_args['projectName']+'_thresholds.pdf')
		pdf_merger.append('Parameter_Definitions.pdf')
		pdf_merger.write(outdir + '/'+pipeline_args['projectName']+'_final_report.pdf')
		pdf_merger.close()

		
		# remove extraneous intermediate files
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_cover_page.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_bulk_data.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_summary_page.pdf')
		stage_for_deletion.append(outdir + '/'+pipeline_args['projectName']+'_thresholds.pdf')

		# actually remove files in stage_for_deletion
		print '\n\n\n' + "Cleaning up project directory"
		for files in stage_for_deletion:
			subprocess.call(['rm', '-rf', files])
