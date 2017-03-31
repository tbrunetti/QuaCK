from fpdf import FPDF
import re
import pandas
import matplotlib.pyplot as plt
import statistics as stats
import collections
import numpy as np

def illumina_snp_overview(inputFile, pdf, clusterSep, aatmean, aatdev, bbtmean, bbtdev, aarmean, abrmean, bbrmean, callrate, outDir):
	print "Running Illumina SNP QC"
	snps_to_remove_text = open(outDir+'/snps_to_remove.txt', 'w')
	snp_qc_table = pandas.read_table(inputFile)
	autosomes_only = snp_qc_table.loc[snp_qc_table['Chr'].isin(list(range(1,23)))] # make a dataframe for chromosomes 1-22 
	non_autosomes = snp_qc_table.loc[~snp_qc_table['Chr'].isin(list(range(0,23)))] # make a dataframe for MT, X, Y, X*Y
	missing_chr = snp_qc_table.loc[snp_qc_table['Chr'] == 0] # make a dataframe for chr 0 only (missing chr values)
	
	
	# calculates total number of snps analzyed and breaks down by category
	total_snps = len(list(snp_qc_table['Name']))
	total_autosomes = len(list(autosomes_only['Name']))
	total_non_autosomes = len(list(non_autosomes['Name']))
	total_missing_chr = len(list(missing_chr['Name']))
	
	columns_for_analysis = []
	print "		Extracting headers needed for analysis"
	# extract all header column names need to SNP analysis	
	for header in list(snp_qc_table):
		if re.search('(.*(Cluster Sep))', header):
			clus_sep = re.search('(.*(Cluster Sep))', header)
			columns_for_analysis.append(clus_sep.group(0))
		
		elif re.search('(.*(AA T Mean))', header):
			AATmean = re.search('(.*(AA T Mean))', header)
			columns_for_analysis.append(AATmean.group(0))

		elif re.search('(.*(AA T Dev))', header):
			AATdev = re.search('(.*(AA T Dev))', header)
			columns_for_analysis.append(AATdev.group(0))

		elif re.search('(.*(BB T Mean))', header):
			BBTmean = re.search('(.*(BB T Mean))', header)
			columns_for_analysis.append(BBTmean.group(0))

		elif re.search('(.*(BB T Dev))', header):
			BBTdev = re.search('(.*(BB T Dev))', header)
			columns_for_analysis.append(BBTdev.group(0))

		elif re.search('(.*(AA R Mean))', header):
			AARmean = re.search('(.*(AA R Mean))', header)
			columns_for_analysis.append(AARmean.group(0))

		elif re.search('(.*(AB R Mean))', header):
			ABRmean = re.search('(.*(AB R Mean))', header)
			columns_for_analysis.append(ABRmean.group(0))

		elif re.search('(.*(BB R Mean))', header):
			BBRmean = re.search('(.*(BB R Mean))', header)
			columns_for_analysis.append(BBRmean.group(0))

	# snp failture stats across ALL chromosomes
	print "		Calculating 5 summary statstics"
	all_stats = {}
	for header_stats in columns_for_analysis:
		all_stats[header_stats] = [stats.median(list(snp_qc_table[header_stats])), stats.mean(list(snp_qc_table[header_stats])), stats.stdev(list(snp_qc_table[header_stats])), 
					min(list(snp_qc_table[header_stats])), max(list(snp_qc_table[header_stats]))]

	
	# snp failure stats call rate, autosomes 
	stats_auto_snp_callrate = [stats.median(list(autosomes_only['Call Freq'])), stats.mean(list(autosomes_only['Call Freq'])), stats.stdev(list(autosomes_only['Call Freq'])), 
					min(list(autosomes_only['Call Freq'])), max(list(autosomes_only['Call Freq']))]

	print "		Extracting failed SNPs"				
	# cluster separation calculations
	snps_fail_clus_sep = list(snp_qc_table[snp_qc_table[clus_sep.group(0)] <= float(clusterSep)]['Name'])
	total_snps_passing_clust = total_snps - len(snps_fail_clus_sep)

	# AA_T mean score calculations
	snps_fail_AATmean = list(snp_qc_table[snp_qc_table[AATmean.group(0)] > float(aatmean)]['Name'])
	total_snps_passing_AATmean = total_snps - len(snps_fail_AATmean)

	# AA_T dev score calculations
	snps_fail_AATdev = list(snp_qc_table[snp_qc_table[AATdev.group(0)] > float(aatdev)]['Name'])
	total_snps_passing_AATdev = total_snps - len(snps_fail_AATdev)

	# BB T mean score calculations
	snps_fail_BBTmean = list(snp_qc_table[snp_qc_table[BBTmean.group(0)] < float(bbtmean)]['Name'])
	total_snps_passing_BBTmean = total_snps - len(snps_fail_BBTmean)

	# BB T dev score calculations
	snps_fail_BBTdev = list(snp_qc_table[snp_qc_table[BBTdev.group(0)] > float(bbtdev)]['Name'])
	total_snps_passing_BBTdev = total_snps - len(snps_fail_BBTdev)

	# AA R mean score calculations
	snps_fail_AARmean = list(snp_qc_table[snp_qc_table[AARmean.group(0)] <= float(aarmean)]['Name'])
	total_snps_passing_AARmean = total_snps - len(snps_fail_AARmean)

	# AB R mean score calculations
	snps_fail_ABRmean = list(snp_qc_table[snp_qc_table[ABRmean.group(0)] <= float(abrmean)]['Name'])
	total_snps_passing_ABRmean = total_snps - len(snps_fail_ABRmean)

	# BB R mean score calculations
	snps_fail_BBRmean = list(snp_qc_table[snp_qc_table[BBRmean.group(0)] <= float(bbrmean)]['Name'])
	total_snps_passing_BBRmean = total_snps - len(snps_fail_BBRmean)

	# remove autosomal SNPs below acceptable call rate
	autosomal_snps_fail_callrate = list(autosomes_only[autosomes_only['Call Freq'] < float(callrate)]['Name'])
	total_snps_passing_callrate = total_autosomes - len(autosomal_snps_fail_callrate)

	# concatenate all failing SNPs
	snps_to_remove = list(set(snps_fail_clus_sep + snps_fail_AATmean + snps_fail_AATdev + snps_fail_BBTmean + snps_fail_BBTdev + snps_fail_AARmean + snps_fail_ABRmean + snps_fail_BBRmean + autosomal_snps_fail_callrate))

	snps_to_remove_text.write('\n'.join(snps_to_remove))
	snps_to_remove_text.flush()

	

	print "		Writing SNP QC to PDF"	
	# create running title for SNP quality
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)
	pdf.cell(0, 30, "SNP Quality Assessment", 0, 1, 'L')
	pdf.line(10, 32, 210, 32)
	pdf.set_font('Arial', '', 12)
	# writes totals into PDF format
	pdf.cell(0, 8, "Total SNPs analyzed:  "+str(total_snps), 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Total autosomal SNPs analyzed:  "+str(total_autosomes), 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Total non-autosomal SNPs analyzed:  "+str(total_non_autosomes), 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Total missing chromosome ID SNPs analyzed:  "+str(total_missing_chr), 0, 1, 'L')
	pdf.cell(0, 8, "Total SNPs passing QC:  "+str(total_snps - len(snps_to_remove)) + ' ' + 
			'('+str(float(total_snps - len(snps_to_remove))/float(total_snps))+'%)' , 0, 1, 'L')
	
	


	# write cluster sep stats 
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 10, "Cluster Separation Statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing cluster separation threshold:  "+str(total_snps_passing_callrate) + '  ' 
			+ '('+str((float(total_snps_passing_callrate)/float(total_autosomes))*100)+'%)', 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Median autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[0]), 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Mean autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[1]), 0, 1, 'L')
	pdf.set_x(20)	
	pdf.cell(0, 8, "Standard deviation autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[2]), 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Minimum autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[3]), 0, 1, 'L')
	pdf.set_x(20)
	pdf.cell(0, 8, "Maximum autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[4]), 0, 1, 'L')
	
	# write call rate score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "call rate score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total autosomal SNPs passing call rate threshold:  "+str(total_autosomes - len(autosomal_snps_fail_callrate)) + '  ' 
			+ '('+str((float(total_autosomes - len(autosomal_snps_fail_callrate))/float(total_autosomes))*100)+'%)', 0, 1, 'L')


	# write AA_T mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "AA T mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AA T mean threshold:  "+str(total_snps_passing_AATmean) + '  ' 
			+ '('+str((float(total_snps_passing_AATmean)/float(total_snps))*100)+'%)', 0, 1, 'L')


	# write AA_T dev score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "AA T dev score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AA T dev threshold:  "+str(total_snps_passing_AATdev) + '  ' 
			+ '('+str((float(total_snps_passing_AATdev)/float(total_snps))*100)+'%)', 0, 1, 'L')


	# write BB_T mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "BB T mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing BB T mean threshold:  "+str(total_snps_passing_BBTmean) + '  ' 
			+ '('+str((float(total_snps_passing_BBTmean)/float(total_snps))*100)+'%)', 0, 1, 'L')

	# write BB_T dev score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "BB T dev score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing BB T dev threshold:  "+str(total_snps_passing_BBTdev) + '  ' 
			+ '('+str((float(total_snps_passing_BBTdev)/float(total_snps))*100)+'%)', 0, 1, 'L')


	# write AA R mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "AA R mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AA R mean threshold:  "+str(total_snps_passing_AARmean) + '  ' 
			+ '('+str((float(total_snps_passing_AARmean)/float(total_snps))*100)+'%)', 0, 1, 'L')


	# write AB R mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "AB R mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AB R mean threshold:  "+str(total_snps_passing_ABRmean) + '  ' 
			+ '('+str((float(total_snps_passing_ABRmean)/float(total_snps))*100)+'%)', 0, 1, 'L')

	# write BB R mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 8, "BB R mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing BB R mean threshold:  "+str(total_snps_passing_BBRmean) + '  ' 
			+ '('+str((float(total_snps_passing_BBRmean)/float(total_snps))*100)+'%)', 0, 1, 'L')

	print "		...Finished writing SNP QC statistics..."	
	return snps_to_remove_text
