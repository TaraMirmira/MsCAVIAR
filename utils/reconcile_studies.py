# Given a list of MsCAVIAR files for a locus, reconcile them so that they
#  contain the same set of SNPs
import argparse
import subprocess
import sys


def reconcile_parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given a list of MsCAVIAR " +
		" LD files for a locus, take out the SNPs that aren't in all studies.")
	parser.add_argument('--infiles', required=True, nargs='+',
		help = 'Space-separated list of MsCAVIAR LD files to LD prune.')
	parser.add_argument('--zscore_file_ext', default = '.zscores',
		help = 'Filename extension for zscores files.')
	parser.add_argument('--ld_file_ext', default = '.ld',
		help = 'Filename extension for LD files.')
	parser.add_argument('--processed_file_ext', default = 'NONE',
		help = 'Use this if you want to include PAINTOR .processed files')
	args = parser.parse_args()
	return args


def get_snps_per_study(z_fnames):
	# snps_by_study is a list where each entry is a dict mapping SNPs to their
	#    line number for the corresponding study in z_fnames
	snps_by_study = []
	for fname in z_fnames:
		snps_by_study.append({})
		with(open(fname, 'r')) as infile:
			line_count = 0
			for line in infile:
				line_count += 1
				snp = line.split('\t')[0]
				snps_by_study[-1][snp] = line_count
	return snps_by_study


# Simple method to take the intersect of two dicts of SNPs
def intersect_another_study(snp_intersect, snps_in_next_study):
	new_intersect = {}
	for snp in snp_intersect:
		if snp in snps_in_next_study:  # SNP is in old intersect and new study
			new_intersect[snp] = True
	return new_intersect


# Return a dict of lines to prune out from the MsCAVIAR files, corresponding to
#    the SNPs that aren't in all files (e.g. not in the snp intersect)
def get_snps_to_prune(snps_in_study, snp_intersect):
	prune_list = {}
	for snp in snps_in_study:
		if snp not in snp_intersect:
			prune_list[snps_in_study[snp]] = True  # the snp LINE, not name
	return prune_list


# Rewrite MsCAVIAR locus file without the pruned SNPs
def rewrite_locus_file(zscore_fname, prune_list, header=False):
	outname = zscore_fname + '_TEMP'
	with(open(zscore_fname, 'r')) as infile:
		with(open(outname, 'w')) as outfile:
			if header:
				outfile.write(infile.readline())
			# track which SNP we're on, write if not in list of pruned SNPs
			snp_num = 0
			for line in infile:
				snp_num += 1
				if snp_num not in prune_list:
					outfile.write(line)
	subprocess.Popen(['mv', outname, zscore_fname]).wait()


# Rewrite MsCAVIAR LD file without the pruned SNPs
def rewrite_ld_file(ld_fname, prune_list):
	outname = ld_fname + '_TEMP'
	with(open(ld_fname, 'r')) as infile:
		with(open(outname, 'w')) as outfile:
			# track which SNP we're on, write if not in list of pruned SNPs
			snp_num = 0
			for line in infile:
				snp_num += 1
				if snp_num not in prune_list:
					splits = line.strip().split(' ')
					# exclude SNP columns in prune list as well
					out_snps = [splits[i] for i in range(len(splits)) if i not in prune_list]
					outfile.write(' '.join(out_snps) + '\n')
	subprocess.Popen(['mv', outname, ld_fname]).wait()


def reconcile_main(args = {}):
	if dir(args) == {}:
		args = ld_prune_parseargs()
	if not args.zscore_file_ext.startswith('.'):
		args.zscore_file_ext = '.' + args.zscore_file_ext
	if not args.ld_file_ext.startswith('.'):
		args.ld_file_ext = '.' + args.ld_file_ext
	if args.processed_file_ext != 'NONE' and not args.processed_file_ext.startswith('.'):
		args.processed_file_ext = '.' + args.processed_file_ext

	ld_fnames = args.infiles
	z_fnames = [i.split(args.ld_file_ext)[0] + args.zscore_file_ext for i in ld_fnames]
	snps_by_study = get_snps_per_study(z_fnames)  # get SNP names in all studies

	# now find the intersect of the SNP names in each study: snps in all studies
	# start with a dict with snps in the first study; values not important
	# then intersect the rest of the studies one-by-one
	snp_intersect = snps_by_study[0]
	for i in range(1, len(snps_by_study)):
		snp_intersect = intersect_another_study(snp_intersect, snps_by_study[i])

	for i in range(len(ld_fnames)):
		prune_list = get_snps_to_prune(snps_by_study[i], snp_intersect)
		rewrite_locus_file(z_fnames[i], prune_list)
		rewrite_ld_file(ld_fnames[i], prune_list)
		if not args.processed_file_ext == 'NONE':
			proc_fname = ld_fnames[i].split(args.ld_file_ext)[0] + args.processed_file_ext
			rewrite_locus_file(proc_fname, prune_list, header=True)


if __name__ == '__main__':
	args = reconcile_parseargs()
	reconcile_main(args)
#
