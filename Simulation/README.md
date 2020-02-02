# MsCaviar Simulation Scripts
For reproducibility purposes, we provide our simulation process:

In order to evaluate the performance of MsCAVIAR as compared with other methods, we performed a simulation study. 
In order to select realistic loci for fine-mapping, we identified regions in a trans-ethnic GWAS of rheumatoid arthritis that contained peak SNPs with p-values of less than 0.0001 and contained ten or more SNPs in a 100kbp region centered around that peak. For each such locus, we used the 1000 Genomes project to generate LD matrices for the SNPs at that locus for both European and East Asian populations. Out of these loci, we selected one region with relatively low LD, where 20% of the SNPs have LD equal to or higher than 0.5, and one region with relatively high LD, where 80% of the SNPs have LD equal to or higher than 0.5. These represent easier and more difficult scenarios, respectively, for fine mapping, since LD makes signals more difficult to distinguish. We pruned groups of SNPs that were in perfect LD in one or more of the populations, leaving one SNP for each. If a group of SNPs were in perfect LD in one population, but not the other, we retained the SNP with the highest Z-score in the other population in order to retain the most signal.

Using these LD matrices, we implanted causal SNPs and simulated their effect sizes. In each simulation,  we implanted either  1, 2, or 3 causal SNPs. Each casual SNP’s true non-centrality parameter Λ was drawn according to N(5.2, 0.125^2). Please refer to our paper for more details:

> Identifying Causal Variants by Fine Mapping Across Multiple Studies  
> Nathan LaPierre, Kodi Taraszka, Helen Huang, Rosemary He, Farhad Hormozdiari, Eleazar Eskin  
> bioRxiv 2020.01.15.908517; doi: https://doi.org/10.1101/2020.01.15.908517

## Sample LD files
In the LD_files/ folder, there are 6 .ld files containing LD matrices that we generated from the 1000 Genomes project. For convenience, we renamed them to {population}\_{amount of LD}.ld, but we label the original file name (chromosome and loci info) in the bracket. In our paper, we reported "Low LD" and "High LD" results, which correspond to EASY and SUPERHARD in the following files.
* _ASN_EASY.ld_ (pruned_locus_chr1_bp_114426001_RA_ASN_2014.ld)
* _ASN_DIFF.ld_ (pruned_locus_chr2_bp_204738919_RA_ASN_2014.ld)
* _ASN_SUPERHARD.ld_ (pruned_locus_chr13_bp_40318819_RA_ASN_2014.ld)
* _EURO_EASY.ld_ (pruned_locus_chr1_bp_114426001_RA_EURO_2014.ld)
* _EURO_DIFF.ld_ (pruned_locus_chr2_bp_204738919_RA_EURO_2014.ld)
* _EURO_SUPERHARD.ld_ (pruned_locus_chr13_bp_40318819_RA_EURO_2014.ld)

The 7th file (_ld_levels_report.txt_) reports the amount of LD for each file above (please use the original file name to search for the corresponding LD level).

## Generating summary statistics for multiple populations
The following Python script takes LD matrices and simulates summary statistics for multiple populations with the **same** sample size (Note: v3 guarantees that non-causal snps would be significant between 5%-50% of the time, and that causal snps would have ld less than 70% among them; v4 guarantees that non-causal snps would be significant between 5%-80% of the time, and impose no requirement on the ld between the causal snps. These requirements control the "difficulty" of the simulated studies):
* _ld_simulate_helen_v3.py_
* _ld_simulate_helen_v4.py_

The following Python script takes LD matrices and simulates summary statistics for multiple populations with **unequal** sample sizes (Note: v5 guarantees that non-causal snps would be significant between 5%-80% of the time, and that causal snps would have ld less than 70% among them; v6 guarantees that non-causal snps would be significant more then 5% of the time (no upper-bound), and impose no requirement on the ld between the causal snps. These requirements control the "difficulty" of the simulated studies):
* _ld_simulate_helen_v5.py_
* _ld_simulate_helen_v6.py_

To run a simulation for the same sample size, run the following line:
```
python3 ld_simulate_helen_v3.py -l1 ${ld_dir}/${pop1}.ld -l2 ${ld_dir}/${pop2}.ld -o $sim_dir -c $this_num_causal -s $num_sim -t $this_tau_2
```
To run a simulation for unequal sample sizes, run the following line:
```
python3 ld_simulate_helen_v5.py -l1 ${ld_dir}/${pop1}.ld -l2 ${ld_dir}/${pop2}.ld -o $sim_dir -c $this_num_causal -s $num_sim -t $this_tau_2 -n $study_size_1','$this_study_size_2
```
where ${ld_dir} is the directory of the LD files; $sim_dir is the directory to output the simulated studies; $this_num_causal is the number of implanted causal variants, $num_sim is the number of simulations the user wish to generate; $this_tau_2 is the heterogeneity (0.5 in the paper) between the studies. Please refer to our paper for more details.

## Running fine-mapping methods on the generated LD and summary statistics
In our paper, CAVIAR, MsCAVIAR, PAINTOR, and SuSiE were compared of their sensitivity and set sizes. Please refer to the following Shell scripts for an idea of how we ran each method (Note: the file paths need to be replaced with your own file paths, and our computer cluster might have a different job handling system from yours). The Shell scripts include some comments to enhance understanding.

* _MsCAVIAR_simulation.sh_
* _MsCAVIAR_unequal.sh_

### Helper scripts
_capture.py_ is a Python script that captures the *sensitivity* (Note: in the preprint manuscript, we called this "accuracy") and *set size* of the each causal set outputted by the methods. Sensitivity is calculated as: length(union(outputted_causal_set, true_causal_set))/length(true_causal_set); whereas set size is calculated as: length(outputted_causal_set).

> For example, if the causal set contains 6 snps: "rs5", "rs6", "rs10", "rs12", "rs20", "rs25", and the true causal snps that were implanted when we generated the summary statistics are "rs6", "rs10", "rs13", then the sensitivity is 0.667, because only "rs6", "rs10" are successfully detected in the causal set, and the set size is 6.

Because PAINTOR does not output the causal set directly, but returns the posterior probability for each snp, the _paintor.R_ script harvests the list of posterior probabilities that PAINTOR returns, and converts them to a causal set by selecting the snps with the highest posterior probabilities that sum to >= 95%, which is the standard we used for CAVIAR and MsCAVIAR.

_susie2.R_ is an R script that runs the R-based package SuSiE (Sum of Single Effects). We used the ```susie_rss()``` function.
