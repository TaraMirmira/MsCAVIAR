## Introduction

In the utils/ folder, you will find several scripts that can be used to pick loci from summary statistics files (based on Z-scores and number of SNPs), generate LD matrices for those SNPs based on 1000 Genomes project data, prune SNPs in perfect LD with one another, and run MsCAVIAR on all those loci. The individual scripts for doing this are available, as is a script that runs the full pipeline (sumstats_mscaviar_pipeline.py).

This pipeline script is ideal for relatively straightforward fine-mapping analyses from summary statistics using MsCAVIAR, e.g. selecting loci based on Z-score levels from multiple summary statistics files with well-defined ethnic groups. It handles all the intermediate file generation and coordinates the SNPs between both studies for steps such as locus creation and LD pruning.

Below, we show how to run an example, and explain each of the scripts and their options.

Please see the data_sources.txt file in the utils/ folder for more information on the sources of the data used here (UK Biobank, Biobank Japan, 1000 Genomes project, PAINTOR).


## Setup

Enter the utils folder and run the provided script to download the 1000 Genomes data, as well as a helper script from the PAINTOR repository* that takes in a list of SNPs and a population name and extracts the 1000 Genomes LD for those SNPs in that population.

```
cd utils/
./download_1000genomes.sh
```

For us, this took about 20 minutes to run, but it can vary depending on the speed of the internet connection.

(*) https://github.com/gkichaev/PAINTOR_V3.0/blob/master/PAINTOR_Utilities/CalcLD_1KG_VCF.py


## Running an example

Download some sample data:

```
wget https://ucla.box.com/shared/static/qcikbiwpq9ix11fb9pkvswy4rt7sbs9e.gz 
tar -xvzf qcikbiwpq9ix11fb9pkvswy4rt7sbs9e.gz
rm qcikbiwpq9ix11fb9pkvswy4rt7sbs9e.gz
```

These are the reformatted versions of UK Biobank and Biobank Japan summary statistics files for Type 2 Diabetes, like the ones run in the paper.

Run the pipeline:

```
python sumstats_mscaviar_pipeline.py --infiles reformatted_sumstats/* --outdir fine_mapping/ --populations EAS EUR --population_sizes 191764 361194 --mscaviar_loc ../MsCAVIAR --exclude_chromosome_six --processes 4
```

Feel free to adjust the number of --processes according to the machine running the program. Depending on what you set that value to, this is likely to take about 15-60 minutes to run, the vast majority of which is LD matrix generation. Afterwards, take a look at the locus directories created inside the fine_mapping/ folder, and MsCAVIAR's 80% confidence sets in the fine_mapping/results/ folder.

This script runs the whole pipeline, from picking loci to generating LD matrices to LD pruning to running MsCAVIAR. By default, only loci with at least 10 SNPs in both studies with Z-scores above 3.9 and at least one Z-score above 5.2 will be selected. To see other defaults and the rest of the options, run

```
python sumstats_mscaviar_pipeline.py -h
```

Below, we review the individual scripts for picking loci, generating LD from 1000 Genomes, and LD pruning. The MsCAVIAR step run by the full pipeline script is simply a wrapper around the main MsCAVIAR package. The arguments for the wrapper script are the ones covered below, except for the aforementioned --processes.


## format_sumstats.py

This script is actually not a part of the main pipeline, but it may be useful in formatting your summary statistics files. We cannot make a script that covers every situation, however. In this case, this script requires your summary statistics file to have the following columns: chromosome, BP/position, SNP name / RSIDs, reference allele, and alternate allele. The columns can be delimited by tab, comma, or any other distinct delimiter, which can be specified with --delimiter.

Additionally, we require there to be a column with Z scores, OR a column with effect sizes (betas) and a standard error column, OR a column with odds ratios and a standard error column. If Z scores are not provided, they will be calculated from the beta/stderr or odds-ratio/stderr.

An additional wrinkle is that the same SNP must have the same name between studies to be treated as the same SNP. To help with this, there is an option "--reformat_snp_ids" that will change the SNP ID to chromosome:position format, e.g. "3:15392837".

For the exact argument names to use, run:

```
python format_sumstats.py -h
```


## pick_loci.py

This script takes as input several summary statistics files in the following format (which is the output format of format_sumstats.py):

```
chr pos rsid A0 A1 Zscore 
chr1 220925144 rs2172699 G A -0.39024
chr1 220925286 rs443189 C T 2.256
```

The script selects loci that can be considered amenable to fine mapping. Namely, loci that have SNPs present in all studies with high Zscores, with a number of other high Zscore SNPs in the area (otherwise, fine mapping would not be necessary). 

The script identifies all candidate "peak" SNPs with an (absolute value) Zscore above a user-set --min_peak_zscore. It then greedily creates loci by starting from the highest Zscore and proceeding to the lowest, centering windows of --window_size base pairs around the peak SNP. If a candidate "peak" SNP is already within the locus of another peak SNP, it does not get its own locus. However, overlapping loci are allowed (non-peak SNPs may be within the window of multiple peak SNPs). A SNP only has to reach the --min_peak_zscore in one study, but this can be changed with --require_peak_all.

Within the window of a peak SNP, all SNPs with (absolute value) Zscores above --min_snp_zscore are included. This threshold should be less than or equal to --min_peak_zscore. These SNPs may not be worth centering a locus around, but their presence and LD structure may help with fine mapping. SNPs are not included unless they are present in both studies and meet the --min_snp_zscore. Loci with fewer than --min_snps number of SNPs that meet this threshold are discarded.

Finally, we allow users to optionally --exclude_chromosome_six, since this has many Human Leukocyte Antigen (HLA) regions, which cause problems for fine mapping.

For the exact argument names to use, run:

```
python pick_loci.py -h
```


## generate_ld_and_mscaviar_files.py

This script takes a formatted summary statistics file (see above) and the continental population label for the people in this study (options: {AFR,AMR,EAS,EUR,SAS}), and generates an LD matrix based on the 1000 Genomes project, with the help of the CalcLD_1KG_VCF.py helper script from the PAINTOR repository. It then creates the properly-formatted Z-scores file for MsCAVIAR.

Users can specify the --chromosome, but this can usually be read easily from the input file. If you followed the setup scripts above, you will not need to use --helper_script or --thousand_genomes, but in principle you could put in custom filepaths. The --output_prefix is used because the script generates both a .ld file and a .processed file, the latter of which is then used to generate a .zscores file for MsCAVIAR.

For the exact argument names to use, run:

```
python generate_ld_and_mscaviar_files.py -h
```


## ld_prune.py

This script takes as input the LD files in a directory and prunes SNPs in perfect LD in both studies. This may eventually be extended to prune SNPs in perfect LD in only one study, though we do not currently do this because the other studies provide different information for these SNPs. The main parameter is --threshold, which allows pruning of SNPs at less than 1.0 LD, e.g. it could be set to 0.99 or 0.98 or something else. However, it's probably best to leave at its default of 1.0. 

Users can also modify the expected file extensions of .zscores and .ld using --zscore_file_ext and --ld_file_ext, but there is generally not much value in this. The different Z scores and LD files must have the same file extension, and the Z score and LD files for a single study should have the same file name (except for the extension). The --processed_file_ext argument can be set to .processed if you also want to LD prune the .processed files generated by generate_ld_and_mscaviar_files.py, which is useful if you want to run PAINTOR later, because PAINTOR uses those files.

For the exact argument names to use, run:

```
python ld_prune.py -h
```

