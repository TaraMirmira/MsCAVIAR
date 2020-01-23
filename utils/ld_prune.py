# Given a list of MsCAVIAR files for a locus, LD prunes SNPs above a threshold
import argparse, sys


def ld_prune_parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given a list of MsCAVIAR " +
		" LD files for a locus, LD prunes SNPs above a threshold.")
	parser.add_argument('--infiles', required=True, nargs='+',
		help = 'Space-separated list of MsCAVIAR LD files to LD prune.')
	parser.add_argument('--output_prefix', required=True,
		help = 'Prefix for pruned output files.')
	parser.add_argument('--threshold', default=1.0, type=float,
		help = 'LD prune SNPs above this threshold of LD with each other.')
	#parser.add_argument('--any_study', action='store_true',
	#       help = 'LD prune SNPs above --threshold in ANY study (default is ALL).')
	#parser.add_argument('--any_pick_criteria', default='max_overall',
	#       choices = ['max_overall', 'max_average', 'max_difference'],
	#       help = 'For any_study, pick SNP with "max" zscore according to this critera.')
	parser.add_argument('--zscore_file_ext', default = '.zscores',
		help = 'Filename extension for zscores files.')
	parser.add_argument('--ld_file_ext', default = '.ld',
		help = 'Filename extension for LD files.')
	parser.add_argument('--processed_file_ext', default = 'NONE',
		help = 'Use this if you want to include PAINTOR .processed files')
	args = parser.parse_args()
	return args


# Find all SNP pairs that are above threshold ("match") in any study
def find_snp_match_by_study(all_infiles, threshold):
	snp_match_by_study = {}  # SNP --> matching SNPs --> studies they match in
	for ld_fname in all_infiles:
		with(open(ld_fname, 'r')) as infile:
			row_num = 0
			for line in infile:
				# symmetric matrix -- only need to check upper triangle
				# store that any SNP pair above threshold matches for this file
				all_corr_for_snp = [float(i) for i in line.strip().split()]
				for col_num in range(row_num + 1, len(all_corr_for_snp)):
					if all_corr_for_snp[col_num] >= threshold:
						if row_num not in snp_match_by_study:
							snp_match_by_study[row_num] = {}
						if col_num not in snp_match_by_study[row_num]:
							snp_match_by_study[row_num][col_num] = []
						snp_match_by_study[row_num][col_num].append(ld_fname)
				row_num += 1
	return snp_match_by_study


# Find "clusters" of SNPs above LD threshold with each other in all studies
def make_clusters_all_study(snp_match_by_study, num_studies):
	clusters = []
	sorted_snp_list = sorted([i for i in snp_match_by_study])
	for snp in sorted_snp_list:
		this_snp_cluster = [snp]
		# find snps that "match" (above LD threshold) in all studies
		for match_snp in snp_match_by_study[snp]:
			if len(snp_match_by_study[snp][match_snp]) == num_studies:
				this_snp_cluster.append(match_snp)
		if len(this_snp_cluster) > 1:  # other SNPs match this one in all studies
			clusters.append(sorted(this_snp_cluster))
	return clusters


'''
# Find "clusters" of SNPs above LD threshold with each other in any study
def make_clusters_any_study(snp_match_by_study):
	clusters = []
	sorted_snp_list = sorted([i for i in snp_match_by_study])
	# store all matching SNP pairs as clusters
	for snp in sorted_snp_list:
		this_snp_cluster = [snp]
		for match_snp in snp_match_by_study[snp]:
			this_snp_cluster.append(match_snp)
		clusters.append(sorted(this_snp_cluster))
	# now merge clusters that overlap, since SNPs in LD with each other
	return clusters
'''


# Pick one SNP to keep from each cluster; mark the rest to be pruned out.
# In this case, we keep the first SNP per cluster, since the choice is arbitrary
def create_prune_list_all_study(clusters):
	prune_list = {}
	for i in range(len(clusters)):
		for j in range(1, len(clusters[i])):
			prune_list[clusters[i][j]] = True
	return prune_list


# Rewrite MsCAVIAR locus file without the pruned SNPs
def rewrite_locus_file(output_prefix, zscore_fname, prune_list, header=False):
	outname = output_prefix + zscore_fname.split('/')[-1]
	with(open(zscore_fname, 'r')) as infile:
		with(open(outname, 'w')) as outfile:
			if header:
				outfile.write(infile.readline())
			# track which SNP we're on, write if not in list of pruned SNPs
			snp_num = 0
			for line in infile:
				if snp_num not in prune_list:
					outfile.write(line)
				snp_num += 1


# Rewrite MsCAVIAR LD file without the pruned SNPs
def rewrite_ld_file(output_prefix, ld_fname, prune_list):
	outname = output_prefix + ld_fname.split('/')[-1]
	with(open(ld_fname, 'r')) as infile:
		with(open(outname, 'w')) as outfile:
			# track which SNP we're on, write if not in list of pruned SNPs
			snp_num = 0
			for line in infile:
				if snp_num not in prune_list:
					splits = line.strip().split(' ')
					# exclude SNP columns in prune list as well
					out_snps = [splits[i] for i in range(len(splits)) if i not in prune_list]
					outfile.write(' '.join(out_snps) + '\n')
				snp_num += 1


def ld_prune_main(args = {}):
	if dir(args) == {}:
		args = ld_prune_parseargs()
	if not args.zscore_file_ext.startswith('.'):
		args.zscore_file_ext = '.' + args.zscore_file_ext
	if not args.ld_file_ext.startswith('.'):
		args.ld_file_ext = '.' + args.ld_file_ext
	if args.processed_file_ext != 'NONE' and not args.processed_file_ext.startswith('.'):
		args.processed_file_ext = '.' + args.processed_file_ext
	snp_match_by_study = find_snp_match_by_study(args.infiles, args.threshold)
	#if args.any_study:
	#	clusters = make_clusters_any_study(snp_match_by_study)
	#	prune_list = create_prune_list_any_study(clusters, args.any_pick_criteria)
	#else:
	num_studies = len(args.infiles)
	clusters = make_clusters_all_study(snp_match_by_study, num_studies)
	prune_list = create_prune_list_all_study(clusters)
	for ld_fname in args.infiles:
		# infer zscore filename from LD filename
		zscore_fname = ld_fname.split(args.ld_file_ext)[0] + args.zscore_file_ext
		rewrite_locus_file(args.output_prefix, zscore_fname, prune_list)
		rewrite_ld_file(args.output_prefix, ld_fname, prune_list)
		if not args.processed_file_ext == 'NONE':
			proc_fname = ld_fname.split(args.ld_file_ext)[0] + args.processed_file_ext
			rewrite_locus_file(args.output_prefix, proc_fname, prune_list, header=True)


if __name__ == '__main__':
	args = ld_prune_parseargs()
	ld_prune_main(args)
#
