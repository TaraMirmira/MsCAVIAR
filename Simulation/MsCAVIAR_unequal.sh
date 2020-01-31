#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_MsCaviar_helenhua
# $ -cwd = run from this current working directory (relative paths)
# -o stdout-Simulation.out
# $ -l h_data=6G,h_rt=24:00:00
# $ -t 1-10:1

# source ~/.bash_profile
# source ~/.bashrc
. /u/local/Modules/default/init/modules.sh
module load python/anaconda3
module load R
# SGE_TASK_ID=1
num_sim=$1
num_repeat=$2
sim_init_dir=$3 # the directory of simulation files
ld_dir=$4 # the directory of ld files
out_init_dir=$5 # output directory
pop1='ASN_EASY'
pop2='EURO_EASY'
# pop1='ASN_DIFF'
# pop2='EURO_DIFF'
# pop1='british10k_1_pruned'
# pop2='asianAll_1_pruned'

# tau_2=(0 0.5 2 0 0.5 2 0 0.5 2)
# num_causal=(1 1 1 2 2 2 3 3 3)
# this_tau_2=${tau_2[$SGE_TASK_ID-1]}
# this_num_causal=${num_causal[$SGE_TASK_ID-1]}
this_tau_2=0.5
this_num_causal=3
study_size_1=10000
study_size_2=(10000 20000 50000 100000)
this_study_size_2=${study_size_2[$SGE_TASK_ID-1]}

if [ ! -d ${sim_init_dir} ]
then
  mkdir ${sim_init_dir}
fi

if [ ! -d ${out_init_dir} ]
then
  mkdir ${out_init_dir}
fi

info='n='$this_study_size_2
echo ${out_init_dir}'_'${info}
mkdir ${sim_init_dir}/${info}
mkdir ${out_init_dir}/${info}

for j in $(seq 1 $num_repeat); do
    sim_dir=${sim_init_dir}/${info}/$j
    out_dir=${out_init_dir}/${info}/$j

    mkdir $sim_dir
    mkdir $out_dir
    
    # if don't want to simulate again, just comment this out
    python3 /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Simulations/ld_simulate_helen_v5.py -l1 ${ld_dir}/${pop1}.ld -l2 ${ld_dir}/${pop2}.ld -o $sim_dir -c $this_num_causal -s $num_sim -t $this_tau_2 -n $study_size_1','$this_study_size_2

    for k in $(seq 1 $num_sim); do
        # for MsCaviar, generate .txt files with paths of ld and zscore
        echo ${ld_dir}/${pop1}.ld > ${sim_dir}/ld_${k}.txt
        echo ${ld_dir}/${pop2}.ld >> ${sim_dir}/ld_${k}.txt
        echo ${sim_dir}/${pop1}_${k}.caviar > ${sim_dir}/gwas_${k}.txt
        echo ${sim_dir}/${pop2}_${k}.caviar >> ${sim_dir}/gwas_${k}.txt

        # run MsCAVIAR and CAVIAR
        /u/home/h/helenhua/project-zarlab/MsCaviar/MsCAVIAR -l ${sim_dir}/ld_${k}.txt -z ${sim_dir}/gwas_${k}.txt -o ${out_dir}/mscaviar_${info}'_'$k -t $this_tau_2 -c $this_num_causal -n $study_size_1','$this_study_size_2
        /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/CAVIAR -l ${ld_dir}/${pop1}.ld -z ${sim_dir}/${pop1}_${k}.caviar -o ${out_dir}/caviar_${pop1}${info}'_'$k -c $this_num_causal
        /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/CAVIAR -l ${ld_dir}/${pop2}.ld -z ${sim_dir}/${pop2}_${k}.caviar -o ${out_dir}/caviar_${pop2}${info}'_'$k -c $this_num_causal

        # for PAINTOR, fake annotate '0's for all locus
        x=`wc -l ${sim_dir}/${pop1}_${k}.caviar | cut -d ' ' -f 1`
        y=$(( x - 1))
        echo "fake_annotate" > ${sim_dir}/paintor_${k}.txt
        for i in $(seq 0 $y); do
          echo "0" >> ${sim_dir}/paintor_${k}.txt
        done
        # for paintor, assemble multiple pops to the same file
        echo 'pos zscore_pop1 pos zscore_pop2' > ${sim_dir}/bothpops_${k}.locus
        paste ${sim_dir}/${pop1}_${k}.caviar ${sim_dir}/${pop2}_${k}.caviar >> ${sim_dir}/bothpops_${k}.locus
        # replace tab with space and reformat files
        sed -i 's/\t/ /g' ${sim_dir}/bothpops_${k}.locus
        cp ${ld_dir}/${pop1}.ld ${sim_dir}/bothpops_${k}.locus.pop1LD
        cp ${ld_dir}/${pop2}.ld ${sim_dir}/bothpops_${k}.locus.pop2LD
        mv ${sim_dir}/paintor_${k}.txt ${sim_dir}/bothpops_${k}.locus.annotations
        echo "bothpops_${k}.locus" > ${sim_dir}/infiles_paintor_${k}.txt
        if [ ! -d ${out_dir}/PAINTOR_${this_num_causal}/ ]
        then
          mkdir ${out_dir}/PAINTOR_${this_num_causal}/
        fi

        # run PAINTOR
        /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/PAINTOR -input ${sim_dir}/infiles_paintor_${k}.txt -in ${sim_dir}/ -out ${out_dir}/PAINTOR_${this_num_causal}/ -Zhead zscore_pop1,zscore_pop2 -LDname pop1LD,pop2LD -annotations fake_annotate -enumerate $this_num_causal

        # analyze accuracy and set size of CAVIAR, MsCAVIAR, and PAINTOR
        mscaviar_result=${out_dir}/'mscaviar_'${info}'_'$k'_set.txt'
        python3 /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/capture.py -s1 $mscaviar_result -t ${sim_dir}/${pop1}.info -p ${out_init_dir}/${info}/'R_mscaviar_'${info} -w mscaviar
        
        caviar_result1=${out_dir}/'caviar_'${pop1}${info}'_'$k'_set'
        caviar_result2=${out_dir}/'caviar_'${pop2}${info}'_'$k'_set'
        python3 /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/capture.py -s1 $caviar_result1 -s2 $caviar_result2 -t ${sim_dir}/${pop1}.info -p ${out_init_dir}/${info}/'R_caviar_'${info} -w caviar
        
        paintor_result=${out_dir}/'paintor_'${info}'_'$k'_set.txt'
        Rscript /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/paintor.R ${out_dir}/PAINTOR_${this_num_causal}/bothpops_${k}.locus.results ${sim_dir}/${pop1}.info $paintor_result ${out_init_dir}/${info}/'R_paintor_'${info}'_recall_rate.txt' ${out_init_dir}/${info}/'R_paintor_'${info}'_config_size.txt'
        
        # for SuSiE, generate output file paths
        susie_subsets_result1=${out_dir}/'susie_'${pop1}${info}'_'$k'_subsets.txt'
        susie_subsets_result2=${out_dir}/'susie_'${pop2}${info}'_'$k'_subsets.txt'
        
        susie_unionset_result1=${out_dir}/'susie_'${pop1}${info}'_'$k'_set.txt'
        susie_unionset_result2=${out_dir}/'susie_'${pop2}${info}'_'$k'_set.txt'
        
        # run SuSiE
        Rscript /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/susie2.R ${sim_dir}/${pop1}_${k}.caviar ${ld_dir}/${pop1}.ld ${sim_dir}/${pop1}.info $this_num_causal $susie_subsets_result1 $susie_unionset_result1 ${out_init_dir}/${info}/'R_susie_pop1_'${info}'_recall_rate.txt' ${out_init_dir}/${info}/'R_susie_pop1_'${info}'_config_size.txt' ${out_init_dir}/${info}/'R_susie_pop1_'${info}'_num_CS.txt'
        Rscript /u/home/h/helenhua/project-zarlab/MSCAVIAR_SIMULATIONS/Automation/susie2.R ${sim_dir}/${pop2}_${k}.caviar ${ld_dir}/${pop2}.ld ${sim_dir}/${pop2}.info $this_num_causal $susie_subsets_result2 $susie_unionset_result2 ${out_init_dir}/${info}/'R_susie_pop2_'${info}'_recall_rate.txt' ${out_init_dir}/${info}/'R_susie_pop2_'${info}'_config_size.txt' ${out_init_dir}/${info}/'R_susie_pop2_'${info}'_num_CS.txt'

    done
done
