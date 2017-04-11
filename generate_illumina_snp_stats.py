from fpdf import FPDF
import re
import pandas
import matplotlib.pyplot as plt
import statistics as stats
import collections
import numpy as np
import seaborn as sns

def illumina_snp_overview(inputFile, pdf, clusterSep, aatmean, aatdev, bbtmean, bbtdev, aarmean, abrmean, bbrmean, callrate, outDir):
	print "Running Illumina SNP QC"
	snps_to_remove_text = open(outDir+'/snps_to_remove.txt', 'w')
	snp_qc_table = pandas.read_table(inputFile)
	snp_qc_table[['Chr']] = snp_qc_table[['Chr']].astype(str)
	chr_strings = [str(chrm) for chrm in range(1, 23)]
	chr_strings_nonauto = [str(chrm) for chrm in range(0, 23)]
	autosomes_only = snp_qc_table.loc[snp_qc_table['Chr'].isin(chr_strings)] # make a dataframe for chromosomes 1-22 
	non_autosomes = snp_qc_table.loc[~snp_qc_table['Chr'].isin(chr_strings_nonauto)] # make a dataframe for MT, X, Y, X*Y
	missing_chr = snp_qc_table.loc[snp_qc_table['Chr'] == '0'] # make a dataframe for chr 0 only (missing chr values)
	chrX_only = non_autosomes.loc[non_autosomes['Chr'].isin(['X'])] #only make dataframe of chrX SNPs
	
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

	
	# TESTING!!!!!!!  CALL RATE FOR X, Y, MITO FILTERING


	non_auto_boxplot = non_autosomes.boxplot(column='Call Freq', by='Chr')
	plt.show()
	
	a = chrX_only.boxplot(column='Call Freq', by='Chr', return_type='both', showfliers=True)
	plt.show()

	b = sns.boxplot(x='Chr', y='Call Freq', data=chrX_only)
	sns.plt.show()

	autosomes_only_boxplot = sns.boxplot(x='Chr', y='Call Freq', data=autosomes_only)
	plt.show()


	# remove non-autosomal SNPs below acceptable call rate
	xChrm_snps_fail_callrate = list(chrX_only[chrX_only['Call Freq'] < float(callrate)]['Name'])
	total_xchrm_snps_passing_callrate = len(list(chrX_only['Name'])) - len(xChrm_snps_fail_callrate)
	print total_xchrm_snps_passing_callrate
	print len(list(chrX_only['Name']))


	# END OF TESTING!!!!!!!





	# snp failture stats across ALL chromosomes
	print "		Calculating five-number summary statstics"
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


	# remove X chr SNPs below acceptable call rate
	xChrm_snps_fail_callrate = list(chrX_only[chrX_only['Call Freq'] < float(callrate)]['Name'])
	total_xchrm_snps_passing_callrate = len(list(chrX_only['Name'])) - len(xChrm_snps_fail_callrate)
	print total_xchrm_snps_passing_callrate
	print len(list(chrX_only['Name']))


	# concatenate all failing SNPs
	snps_to_remove = list(set(snps_fail_clus_sep + snps_fail_AATmean + snps_fail_AATdev + snps_fail_BBTmean + snps_fail_BBTdev + snps_fail_AARmean + snps_fail_ABRmean + snps_fail_BBRmean + autosomal_snps_fail_callrate + xChrm_snps_fail_callrate))

	snps_to_remove_text.write('\n'.join(snps_to_remove))
	snps_to_remove_text.flush()

	
	# calculate passing SNP stats for top table on page
	autosomes_remain_table = autosomes_only.loc[~autosomes_only['Name'].isin(snps_to_remove)]
	autosomes_remain = len(autosomes_remain_table.index)

	non_auto_remain_table = non_autosomes.loc[~non_autosomes['Name'].isin(snps_to_remove)]
	non_auto_remain = len(non_auto_remain_table.index)

	missing_remain_table = missing_chr.loc[~missing_chr['Name'].isin(snps_to_remove)]
	missing_remain = len(missing_remain_table.index)




	print "		Writing SNP QC to PDF"	
	# create running title for SNP quality
	pdf.add_page()
	pdf.set_font('Arial', 'B', 30)
	pdf.cell(0, 30, "SNP Quality Assessment", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_fill_color(200)
	pdf.set_font('Arial', 'B', 14)
	# writes totals into PDF format
	pdf.multi_cell(0, 8, "Total SNPs analyzed:  "+str(total_snps), 1, 'L', True)
	pdf.set_font('Arial', '', 14)
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total autosomal SNPs analyzed:  "+str(total_autosomes), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total non-autosomal SNPs analyzed:  "+str(total_non_autosomes), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Total missing chromosome ID SNPs analyzed:  "+str(total_missing_chr), 1, 1, 'L')
	pdf.set_font('Arial', 'B', 14)
	pdf.multi_cell(0, 8, "Total SNPs passing QC:  "+str(total_snps - len(snps_to_remove)) + ' ' + 
			'('+str((float(total_snps - len(snps_to_remove))/float(total_snps)) * 100)+'%)' , 1, 'L', True)
	pdf.set_font('Arial', '', 14)
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Autosomal SNPs remaining:  "+str(autosomes_remain) + ' (' +str((float(autosomes_remain)/float(total_autosomes))*100)+'%)', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Non-autosomal SNPs remaining:  "+str(non_auto_remain) + ' (' +str((float(non_auto_remain)/float(total_non_autosomes))*100)+'%)', 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 8, "Missing chromosome ID SNPs remaining:  "+str(missing_remain) + ' (' +str((float(missing_remain)/float(total_missing_chr))*100)+'%)', 1, 1, 'L')
	pdf.multi_cell(0, 8, '\n\n\n', 0, 1, 'L')
	


	# write cluster sep stats 
	pdf.set_font('Arial', 'B', 14)
	pdf.multi_cell(0, 10, "Cluster Separation Statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.multi_cell(0, 8, "Total SNPs passing cluster separation threshold:  "+str(total_snps_passing_clust) + '  ' 
			+ '('+str((float(total_snps_passing_clust)/float(total_snps))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median cluster separation:  "+ str(all_stats[clus_sep.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean cluster separation:  "+ str(all_stats[clus_sep.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of cluster separation:  "+ str(all_stats[clus_sep.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum cluster separation:  "+ str(all_stats[clus_sep.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum cluster sep:  "+ str(all_stats[clus_sep.group(0)][4]), 0, 1, 'L')
	
	
	# write call rate score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "Call rate score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total autosomal SNPs passing call rate threshold:  "+str(total_autosomes - len(autosomal_snps_fail_callrate)) + '  ' 
			+ '('+str((float(total_autosomes - len(autosomal_snps_fail_callrate))/float(total_autosomes))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum autosomal SNP call rate:  "+ str(stats_auto_snp_callrate[4]), 0, 1, 'L')

	
	# call rate statistics for X chromosome
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "Chromosome X call rate score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total chromosome X SNPs passing call rate threshold:  "+str(total_xchrm_snps_passing_callrate) + '  ' 
			+ '('+str((float(total_xchrm_snps_passing_callrate)/float(len(list(chrX_only['Name']))))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median chrX SNP call rate:  "+ str(stats.median(list(chrX_only['Call Freq']))), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean chrX SNP call rate:  "+ str(stats.mean(list(chrX_only['Call Freq']))), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation autosomal chrX call rate:  "+ str(stats.stdev(list(chrX_only['Call Freq']))), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum chrX SNP call rate:  "+ str(min(list(chrX_only['Call Freq']))), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum chrX SNP call rate:  "+ str(max(list(chrX_only['Call Freq']))), 0, 1, 'L')


	

	# write AA_T mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "AA T mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AA T mean threshold:  "+str(total_snps_passing_AATmean) + '  ' 
			+ '('+str((float(total_snps_passing_AATmean)/float(total_snps))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized AA theta mean:  "+ str(all_stats[AATmean.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized AA theta mean:  "+ str(all_stats[AATmean.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized AA theta mean:  "+ str(all_stats[AATmean.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized AA theta mean:  "+ str(all_stats[AATmean.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized AA theta mean:  "+ str(all_stats[AATmean.group(0)][4]), 0, 1, 'L')


	# write AA_T dev score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "AA T dev score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AA T dev threshold:  "+str(total_snps_passing_AATdev) + '  ' 
			+ '('+str((float(total_snps_passing_AATdev)/float(total_snps))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized AA theta deviation:  "+ str(all_stats[AATdev.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized AA theta deviation:  "+ str(all_stats[AATdev.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized AA theta deviation:  "+ str(all_stats[AATdev.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized AA theta deviation:  "+ str(all_stats[AATdev.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized AA theta deviation:  "+ str(all_stats[AATdev.group(0)][4]), 0, 1, 'L')



	# write BB_T mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "BB T mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing BB T mean threshold:  "+str(total_snps_passing_BBTmean) + '  ' 
			+ '('+str((float(total_snps_passing_BBTmean)/float(total_snps))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized BB theta mean:  "+ str(all_stats[BBTmean.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized BB theta mean:  "+ str(all_stats[BBTmean.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized BB theta mean:  "+ str(all_stats[BBTmean.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized BB theta mean:  "+ str(all_stats[BBTmean.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized BB theta mean:  "+ str(all_stats[BBTmean.group(0)][4]), 0, 1, 'L')



	# write BB_T dev score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "BB T dev score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing BB T dev threshold:  "+str(total_snps_passing_BBTdev) + '  ' 
			+ '('+str((float(total_snps_passing_BBTdev)/float(total_snps))*100)+'%)', 0, 1, 'L')

	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized BB theta deviation:  "+ str(all_stats[BBTdev.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized BB theta deviation:  "+ str(all_stats[BBTdev.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized BB theta deviation:  "+ str(all_stats[BBTdev.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized BB theta deviation:  "+ str(all_stats[BBTdev.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized BB theta deviation:  "+ str(all_stats[BBTdev.group(0)][4]), 0, 1, 'L')



	# write AA R mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "AA R mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AA R mean threshold:  "+str(total_snps_passing_AARmean) + '  ' 
			+ '('+str((float(total_snps_passing_AARmean)/float(total_snps))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized AA intesity mean:  "+ str(all_stats[AARmean.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized AA intesity mean:  "+ str(all_stats[AARmean.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized AA intesity mean:  "+ str(all_stats[AARmean.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized AA intesity mean:  "+ str(all_stats[AARmean.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized AA intesity mean:  "+ str(all_stats[AARmean.group(0)][4]), 0, 1, 'L')




	# write AB R mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "AB R mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing AB R mean threshold:  "+str(total_snps_passing_ABRmean) + '  ' 
			+ '('+str((float(total_snps_passing_ABRmean)/float(total_snps))*100)+'%)', 0, 1, 'L')
	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized AB intesity mean:  "+ str(all_stats[ABRmean.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized AB intesity mean:  "+ str(all_stats[ABRmean.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized AB intesity mean:  "+ str(all_stats[ABRmean.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized AB intesity mean:  "+ str(all_stats[ABRmean.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized AB intesity mean:  "+ str(all_stats[ABRmean.group(0)][4]), 0, 1, 'L')


	# write BB R mean score stats
	pdf.set_font('Arial', 'B', 14)
	pdf.cell(0, 15, "BB R mean score statistics", 0, 1, 'L')
	pdf.set_font('Arial', '', 12)
	pdf.cell(0, 8, "Total SNPs passing BB R mean threshold:  "+str(total_snps_passing_BBRmean) + '  ' 
			+ '('+str((float(total_snps_passing_BBRmean)/float(total_snps))*100)+'%)', 0, 1, 'L')

	pdf.multi_cell(0, 8, "Summary Stats on Original Data:")
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Median normalized BB intesity mean:  "+ str(all_stats[BBRmean.group(0)][0]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Mean normalized BB intesity mean:  "+ str(all_stats[BBRmean.group(0)][1]), 0, 1, 'L')
	pdf.set_x(40)	
	pdf.multi_cell(0, 5, "Standard deviation of normalized BB intesity mean:  "+ str(all_stats[BBRmean.group(0)][2]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Minimum normalized BB intesity mean:  "+ str(all_stats[BBRmean.group(0)][3]), 0, 1, 'L')
	pdf.set_x(40)
	pdf.multi_cell(0, 5, "Maximum normalized BB intesity mean:  "+ str(all_stats[BBRmean.group(0)][4]), 0, 1, 'L')


	print "		...Finished writing SNP QC statistics..."	
	return snps_to_remove_text
