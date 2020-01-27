# Given intersected sumstats files, picks loci for fine mapping based on user
#   specified parameters such as min. zscore and number of SNPs.
import argparse, os, sys


def pick_loci_parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given sumstats files, picks" +
		" loci for fine mapping based on zscores and number of SNPs.")
	parser.add_argument('--infiles', required=True, nargs='+',
		help = 'Space-separated list of formatted sumstats files.')
	parser.add_argument('--outdir', required=True,
		help = 'Directory to write locus files in.')
	parser.add_argument('--exclude_chromosome_six', action='store_true',
		help = 'Exclude chromosome six from analysis due to HLA regions.')
	parser.add_argument('--include_non_autosomal', action='store_true',
		help = 'Include non autosomal (X, Y, MT) chromosomes.')
	parser.add_argument('--min_peak_zscore', default=5.2, type=float,
		help = 'Minimum (abs) zscore of peak SNP to create a locus.')
	parser.add_argument('--min_snp_zscore', default=3.9, type=float,
		help = 'Minimum (abs) zscore of SNP to include in a locus.')
	parser.add_argument('--min_snps', default=10, type=int,
		help = 'Minimum number of SNPs in a locus to use it.')
	parser.add_argument('--require_peak_all', action='store_true',
		help = 'Require a peak SNP to be present in a locus in ALL studies.')
	parser.add_argument('--window_size', default=100000, type=int,
		help = 'Window size to center around peak SNP to create a locus.')
	args = parser.parse_args()
	return args


# Find all SNPs above the min_peak_zscore threshold
def find_ordered_peak_snps(all_infiles, min_peak_zscore):
	peak_snps, duplicate_snps = [], {}
	for fname in all_infiles:
		snps_seen = {}
		with(open(fname, 'r')) as infile:
			header = infile.readline()
			for line in infile:
				splits = line.strip().split()
				snp_name, zscore = splits[2], abs(float(splits[-1]))
				if snp_name in snps_seen:
					duplicate_snps[snp_name] = True
				else:
					snps_seen[snp_name] = True
				if zscore >= min_peak_zscore:
					splits[-1] = zscore
					peak_snps.append(splits)
	return peak_snps, duplicate_snps


# Greedily create loci in order of highest abs(zscore), with a window of
#    window_size BP centered around peak. Overlapping loci are kept separate.
#    Peak SNPs in an already-existing locus do not form a new locus.
def pick_locus_locations_greedily(window_size, exclude6, include_sex, peak_snps):
	loci = {}  # dict so we can search by chromosome
	peak_snps.sort(key=lambda x: x[-1], reverse=True)  # sort by zscore
	for snp in peak_snps:
		chrom, pos = snp[0], int(snp[1])
		clean_chrom = chrom.strip('chromosomeCHROMOSOME')
		if exclude6 and clean_chrom == '6':
			continue
		if not include_sex and clean_chrom.upper() in ['X', 'Y', 'MT']:
			continue
		half_window = int(window_size / 2.0)
		start, end = pos - half_window, pos + half_window

		if chrom not in loci:
			loci[chrom] = [[start, end, []]]  # third element: list of SNPs in locus
		else:  # if chromosome already has loci, check if SNP is already in a locus
			inside_another_locus = False
			for window in loci[chrom]:
				if pos > window[0] and pos < window[1]:
					inside_another_locus = True
					break
			if not inside_another_locus:
				loci[chrom].append([start, end, []])

	for chrom in loci:  # now sort loci by position within each chromosome
		loci[chrom].sort(key=lambda x: x[0])
	return loci


# Actually gathers the SNP information for all SNPs in all studies that fall
#    within the locus locations.
def gather_locus_snps(all_infiles, min_snp_zscore, loci, duplicate_snps):
	# loci will now store a 4d array per locus per chromosome:
	#    locus location > study > SNP > info for SNP
	for fname in all_infiles:
		# initialize a list to hold SNPs for this study at each locus
		# if this list remains empty, the study has no SNPs at that locus
		for chrom in loci:
			for i in range(len(loci[chrom])):
				loci[chrom][i][-1].append([])

		with(open(fname, 'r')) as infile:
			header = infile.readline()
			for line in infile:
				splits = line.strip().split()
				chrom, pos, rsid, a0, a1, zscore = splits
				if chrom not in loci:
					continue
				if abs(float(zscore)) < min_snp_zscore:
					continue
				if rsid in duplicate_snps:
					continue

				intpos = int(pos)
				# since loci are ordered, skip SNP if it's before the start of
				#    the first locus and break if it's after end of last locus
				#if intpos < loci[chrom][0][0]:
				#	continue
				#if intpos > loci[chrom][-1][1]:
				#	break
				for i in range(len(loci[chrom])):
					start, end = loci[chrom][i][:2]
					if intpos < start or intpos > end:
						continue
					# else, SNP is in the locus and should be appended to the
					#    list created for this study for this locus (index -1)
					loci[chrom][i][2][-1].append(splits)
	return loci


# Find all SNPs in the first study and establish A0 and A1
def initialize_snp_intersect(study_locus):
	snp_intersect = {}
	for snp_info in study_locus:
		chrom, pos, rsid, a0, a1, zscore = snp_info
		if a0 == a1:
			continue  # exclude monoallelic SNPs
		snp_intersect[rsid] = [a0, a1]
	return snp_intersect


# Intersect another study with the SNPs in the intersect of all previous studies
def intersect_another_study(snp_intersect, study_locus):
	new_intersect = {}
	for snp_info in study_locus:
		chrom, pos, rsid, a0, a1, zscore = snp_info
		if a0 == a1:
			continue  # exclude monoallelic SNPs
		if rsid in snp_intersect:  # now we know it's in all studies so far
			# check to make sure not triallelic or quad-allelic
			if a0 in snp_intersect[rsid] and a1 in snp_intersect[rsid]:
				new_intersect[rsid] = snp_intersect[rsid]
	return new_intersect


# remove all snps from all loci that aren't in the intersect SNPs, and if A0/A1
#    is flipped in a study, flip sign of zscore
def update_locus(locus, snp_intersect):
	new_locus = []
	for study_locus in locus:
		new_study_locus = []
		for snp_info in study_locus:
			chrom, pos, rsid, a0, a1, zscore = snp_info
			if rsid not in snp_intersect:
				continue

			if a0 == snp_intersect[rsid][1]:  # A0/A1 flipped in this study
				zscore = str(-1 * float(zscore))
				a0, a1 = snp_intersect[rsid]
			new_study_locus.append([chrom, pos, rsid, a0, a1, zscore])
		new_locus.append(new_study_locus)
	return new_locus


# Filter out loci that no longer have a suitable peak SNP after intersection.
# False means do not filter out the locus, True means filter out the locus.
def filter_locus(locus, min_peak_zscore, require_peak_all):
	for study in locus:
		study_zscores = [abs(float(snp[-1])) for snp in study]
		max_zscore = max(study_zscores)
		# if a SNP is above the min_peak_zscore and we don't require all studies
		#    to hit the peak, we know this locus doesn't need to be filtered
		if max_zscore >= min_peak_zscore:
			if not require_peak_all:
				return False
		# if we require all studies to have a SNP above min_peak_zscore and this
		#    one doesn't reach the threshold, this locus must be filtered
		elif require_peak_all:
			return True
	# if we never reached the False/True conditions above, return the opposite
	if require_peak_all:
		return False
	else:
		return True


# Write out locus files to the outdir
def write_locus(infiles, outdir, chrom, locus_info):
	start, end, locus = locus_info
	# each locus gets its own directory, in which the locus files for each study
	#    are stored
	locus_name = chrom + '_bp_' + str(start) + '_to_' + str(end)
	locus_outdir = outdir + 'locus_chromosome_' + locus_name + '/'
	if not os.path.exists(locus_outdir):
		os.makedirs(locus_outdir)
	for i in range(len(infiles)):
		outname = locus_outdir + 'study_' + infiles[i].split('/')[-1]
		with(open(outname, 'w')) as outfile:
			outfile.write('chr pos rsid A0 A1 Zscore\n')  # outfile header
			file_locus = locus[i]  # SNPs at the locus for this file
			for snp_info in file_locus:
				outfile.write(' '.join([str(i) for i in snp_info]) + '\n')


def pick_loci_main(args = {}):
	if dir(args) == {}:
		args = pick_loci_parseargs()
	if not args.outdir.endswith('/'):
		args.outdir += '/'
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	# the "min" zscores are really minimum absolute values
	args.min_peak_zscore = abs(args.min_peak_zscore)
	args.min_snp_zscore = abs(args.min_snp_zscore)

	# pick loci and gather the snps at those loci
	peak_snps, duplicate_snps = find_ordered_peak_snps(args.infiles, args.min_peak_zscore)
	loci = pick_locus_locations_greedily(args.window_size,
		args.exclude_chromosome_six, args.include_non_autosomal, peak_snps)
	loci = gather_locus_snps(args.infiles, args.min_snp_zscore, loci, duplicate_snps)

	# intersect SNPs across studies, then filter loci below min_snps or whose
	#    suitable peak SNPs were removed in the intersect step. if not filtered,
	#    write the locus files.
	for chrom in loci:
		for i in range(len(loci[chrom])):
			locus = loci[chrom][i][2]
			snp_intersect = initialize_snp_intersect(locus[0])
			for j in range(1, len(locus)):
				snp_intersect = intersect_another_study(snp_intersect, locus[j])
			if len(snp_intersect) < args.min_snps:
				continue  # not enough SNPs in the SNP intersect
			locus = update_locus(locus, snp_intersect)
			filter_this = filter_locus(locus, args.min_peak_zscore, args.require_peak_all)
			if not filter_this:
				loci[chrom][i][2] = locus
				write_locus(args.infiles, args.outdir, chrom, loci[chrom][i])


if __name__ == '__main__':
	args = pick_loci_parseargs()
	pick_loci_main(args)
#
