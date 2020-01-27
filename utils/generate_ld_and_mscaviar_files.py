# Given a locus file and 1000 genomes location, generates MsCAVIAR format
#    locus file and LD matrix
import argparse, glob, os, subprocess
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))) + '/'


def gen_ld_parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given a locus file and " +
		"1000 genomes location, generates MsCAVIAR format locus file and LD matrix.")
	parser.add_argument('--infile', required=True,
		help = 'Input locus file in format given by format_sumstats.')
	parser.add_argument('--output_prefix', required=True,
		help = 'Output file name prefix.')
	parser.add_argument('--population', required=True,
		choices = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
		help = 'Continental population group for this study.')
	parser.add_argument('--chromosome', default='AUTO',
		help = 'Optionally inform the script which chromosome this is.')
	parser.add_argument('--helper_script', default='AUTO',
		help = 'Helper script from PAINTOR to get 1000genomes LD using infile.')
	parser.add_argument('--thousand_genomes', default='AUTO',
		help = '1000 genomes data location.')
	args = parser.parse_args()
	return args


# attempt to automatically determine locus chromosome by looking at infile
def auto_determine_chromosome(infile):
	with(open(infile, 'r')) as inf:
		header = inf.readline()
		chrom = inf.readline().split()[0]
		return chrom.strip('chromosomeCHROMOSOME')  # try to get just the number


# reformat locus file into MsCAVIAR locus file format
def gen_mscaviar_locus_file(output_prefix):
	with(open(output_prefix + '.processed', 'r')) as infile:
		infile.readline()  # skip header
		with(open(output_prefix + '.zscores', 'w')) as outfile:
			for line in infile:
				splits = line.strip().split()
				rsid, zscore = splits[2], splits[5]
				outfile.write(rsid + '\t' + zscore + '\n')


def gen_ld_main(args = {}):
	if dir(args) == {}:
		args = gen_ld_parseargs()
	# set default locations
	if args.helper_script == 'AUTO':
		args.helper_script = __location__ + '1000genomes/CalcLD_1KG_VCF.py'
	if args.thousand_genomes == 'AUTO':
		args.thousand_genomes = __location__ + '1000genomes/'
	if not args.thousand_genomes.endswith('/'):
		args.thousand_genomes += '/'
	if args.chromosome == 'AUTO':
		args.chromosome = auto_determine_chromosome(args.infile)

	# now run the helper script to get the 1000genomes LD file
	ref_loc = glob.glob(args.thousand_genomes + 'ALL.chr' + args.chromosome + '.*.vcf.gz')[0]
	map_loc = args.thousand_genomes + 'integrated_call_samples_v3.20130502.ALL.panel'
	subprocess.Popen(['python', args.helper_script, '--locus', args.infile,
		'--reference', ref_loc, '--map', map_loc, '--effect_allele', 'A1',
		'--alt_allele', 'A0', '--population', args.population, '--Zhead', 'Zscore',
		'--out_name', args.output_prefix, '--position', 'pos']).wait()
	gen_mscaviar_locus_file(args.output_prefix)


if __name__ == '__main__':
	args = gen_ld_parseargs()
	gen_ld_main(args)
#
