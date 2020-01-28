'''
This is a wrapper script that combines several functions:
- Picking loci from sumstats files based on Zscores and number of SNPs
- Generating LD matrices for these loci based on 1000 Genomes
- LD pruning those loci and intersecting the SNPs between the studies
- Running MsCAVIAR on all those loci
Essentially, the input is multiple sumstats files, ideally from different
    populations, and the output is the MsCAVIAR results for all loci that were
	selected according to user specifications.
'''
import sys, traceback


import argparse, glob, multiprocessing, os, subprocess
import pick_loci as locus_picker
import generate_ld_and_mscaviar_files as ld_generator
import ld_prune as pruner
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))) + '/'


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given formatted SumStats " +
		"files, picks loci, generates LD and MsCAVIAR files, & runs MsCAVIAR.")

	# REQUIRED ARGUMENTS (Sumstats files, population info, MsCAVIAR location)

	parser.add_argument('--infiles', required=True, nargs='+',
		help = 'Space-separated list of formatted sumstats files.')
	parser.add_argument('--outdir', required=True,
		help = 'Directory to write results files to.')
	parser.add_argument('--populations', required=True, nargs = '+',
		choices = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
		help = 'Continental population group for each study.')
	parser.add_argument('--population_sizes', required=True, nargs = '+',
		type = int, help = 'Population size of each study, in order.')

	# LOCUS PICKING ARGUMENTS

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

	# LD GENERATION AND PRUNING ARGUMENTS

	parser.add_argument('--processes', type = int, default = 4,
		help = 'Number of simultaneous processes to run for LD generation.')
	parser.add_argument('--helper_script', default='AUTO',
		help = 'Helper script from PAINTOR to get 1000genomes LD using infile.')
	parser.add_argument('--thousand_genomes', default='AUTO',
		help = '1000 genomes data location.')
	parser.add_argument('--threshold', default=1.0, type=float,
		help = 'LD prune SNPs above this threshold of LD with each other.')

	# MsCAVIAR ARGUMENTS

	parser.add_argument('--mscaviar_loc', default = 'AUTO',
		help = 'MsCAVIAR executable location')
	parser.add_argument('--max_causal', default=3, type=int,
		help = 'Maximum number of causal SNPs allowed per locus')
	parser.add_argument('--rho', default=0.8, type=float,
		help = 'Posterior probability threshold ("rho")')
	parser.add_argument('--tau_sqr', default=0.52, type=float,
		help = 'Heterogeneity parameter ("tau^2")')
	args = parser.parse_args()
	return args


# Figure out which of the input populations this file name belongs to, using the
#    fact that args.infiles is in the same order as args.populations and the
#    base name of an infile will be in fname based on how pick_loci is run.
def determine_population(infiles, fname, populations):
	for i in range(len(infiles)):
		base_name = infiles[i].split('/')[-1].split('.')[0]
		if base_name in fname:
			return populations[i]


def run_ld_generation_process(args, fname):
	try:
		this_args = args
		this_args.infile = fname
		this_args.output_prefix = fname.split('.')[0]
		this_args.population = determine_population(args.infiles, fname, args.populations)
		ld_generator.gen_ld_main(this_args)
		return 0
	except:
		print('Error in process for ' + str(fname))
		err = sys.exc_info()
		tb = traceback.format_exception(err[0], err[1], err[2])
		print(''.join(tb) + '\n')
		sys.stdout.flush()


# For a directory, glob the infiles for the LD pruning script, output to pruned/.
#    Then generate MsCAVIAR -l and -z files and run MsCAVIAR, making sure the
#    original order of studies is maintained.
def run_ld_prune_and_mscaviar_on_dir(args, ldir, mscaviar_results_dir):
	# set up input and output file names
	locus_dir_path = args.outdir + ldir
	locus_fnames = glob.glob(locus_dir_path + '/*.ld')
	this_args = args
	this_args.infiles = locus_fnames
	this_outdir = locus_dir_path + '/pruned/'
	if not os.path.isdir(this_outdir):
		os.makedirs(this_outdir)
	this_args.output_prefix = this_outdir
	# run LD pruning
	pruner.ld_prune_main(this_args)

	# for MsCAVIAR, we must make sure the files are in their original order
	sumstats_fnames = [i.split('/')[-1].split('.')[0] for i in args.infiles]
	ordered_znames = []
	for fname in sumstats_fnames:
		ordered_znames.extend(glob.glob(this_outdir + '*' + fname + '*.zscores'))
	ordered_ldnames = [i.split('.zscores')[0] + '.ld' for i in ordered_znames]
	ordered_popsizes = ','.join([str(i) for i in args.population_sizes])  # MsCAVIAR input format
	# now create the -l and -z files for MsCAVIAR
	zfiles_list_file = this_outdir + 'zfiles.txt'
	ldfiles_list_file = this_outdir + 'ldfiles.txt'
	with(open(zfiles_list_file, 'w')) as outfile:
		for zname in ordered_znames:
			outfile.write(zname + '\n')
	with(open(ldfiles_list_file, 'w')) as outfile:
		for ldname in ordered_ldnames:
			outfile.write(ldname + '\n')
	# finally, run MsCAVIAR
	mscaviar_output_prefix = this_outdir + 'mscaviar_results_' + ldir
	subprocess.Popen([args.mscaviar_loc, '-l', ldfiles_list_file, '-z',
		zfiles_list_file, '-r', str(args.rho), '-t', str(args.tau_sqr),
		'-c', str(args.max_causal), '-n', ordered_popsizes,
		'-o', mscaviar_output_prefix]).wait()
	# copy results into easier-to-find directory
	subprocess.Popen(['cp', mscaviar_output_prefix + '_set.txt',
		mscaviar_results_dir]).wait()


def main():
	args = parseargs()
	if not args.outdir.endswith('/'):
		args.outdir += '/'
	if args.mscaviar_loc == 'AUTO':
		args.mscaviar_loc = __location__ + '../MsCAVIAR'
	args.zscore_file_ext = '.zscores'
	args.ld_file_ext = '.ld'
	args.processed_file_ext = '.processed'
	args.chromosome = 'AUTO'

	# Run pick loci. glob to get all resulting locus files for all studies.
	locus_picker.pick_loci_main(args)
	all_fnames = glob.glob(args.outdir + '*/study_*')

	# Use the above list to run gen_ld script. This is the slowest step but is
	#    also embarassingly parallel, so we use multiprocessing to speed it up.
	if args.processes > 1:
		pool = multiprocessing.Pool(processes = args.processes)
		for fname in all_fnames:
			pool.apply_async(run_ld_generation_process, (args, fname,))
		pool.close()
		pool.join()
	else:
		for fname in all_fnames:
			run_ld_generation_process(args, fname)

	# Get list of locus dirs, run LD pruning and MsCAVIAR on each.
	mscaviar_results_dir = args.outdir + 'results/'
	if not os.path.isdir(mscaviar_results_dir):
		os.makedirs(mscaviar_results_dir)
	locus_dirs = next(os.walk(args.outdir))[1]
	locus_dirs = [i for i in locus_dirs if 'locus_chromosome_' in i]
	for ldir in locus_dirs:
		run_ld_prune_and_mscaviar_on_dir(args, ldir, mscaviar_results_dir)


if __name__ == '__main__':
	main()
#
