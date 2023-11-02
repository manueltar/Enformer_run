#!/bin/bash

analysis=$1
MASTER_ROUTE=$2

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

#rm -rf $output_dir
#mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files


#### export_vcf #############################


type=$(echo "export_vcf""_""$analysis")
outfile_export_vcf=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_export_vcf
echo -n "" > $outfile_export_vcf
name_export_vcf=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_export_vcf=$(echo "$Rscripts_path""31_vcf_exporter.R")


Table_S6=$(echo "/group/soranzo/manuel.tardaguila/Motif_analysis/""Table_S6.tsv")

myjobid_export_vcf=$(sbatch --job-name=$name_export_vcf --output=$outfile_export_vcf --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_export_vcf --Table_S6 $Table_S6 --type $type --out $output_dir")
myjobid_seff_export_vcf=$(sbatch --dependency=afterany:$myjobid_export_vcf --open-mode=append --output=$outfile_export_vcf --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_export_vcf >> $outfile_export_vcf")



#### conda environment for enforcer  #######


 module load cuda/11.5.1
 module load cudnn/8.3.3.40-11.5-cuda-11.5.1
 
 conda activate enformer_20231026



type=$(echo "Enformer_run""_""$analysis")
outfile_Enformer_run=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_Enformer_run
echo -n "" > $outfile_Enformer_run
name_Enformer_run=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


input_vcf=$(echo "$output_dir""VARS.vcf")
variants_list=$(echo "$output_dir""variants_list.txt")
output_folder=$(echo "$output_dir""Enformer_results""/")
#reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2022.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")
reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2023.1/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa")
enformer_predict_script=$(echo "/group/soranzo/manuel.tardaguila/Enformer/enformer_variant_score/enformer_predict.py")



myjobid_Enformer_run=$(sbatch --dependency=afterany:$myjobid_export_vcf --job-name=$name_Enformer_run --output=$outfile_Enformer_run --partition=gpuq --time=24:00:00 --nodes=1 --ntasks=4 --gres=gpu:2 --mem=16GB --parsable --wrap="$enformer_predict_script --vcf $input_vcf --output $output_folder --ref_genome $reference_genome --add_chr_prefix")
myjobid_seff_Enformer_run=$(sbatch --dependency=afterany:$myjobid_Enformer_run --open-mode=append --output=$outfile_Enformer_run --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Enformer_run >> $outfile_Enformer_run")

# myjobid_Enformer_run=$(sbatch --job-name=$name_Enformer_run --output=$outfile_Enformer_run --partition=gpuq --time=24:00:00 --nodes=1 --ntasks=4 --gres=gpu:2 --mem=16GB --parsable --wrap="$enformer_predict_script --vcf $input_vcf --output $output_folder --ref_genome $reference_genome --add_chr_prefix")
# myjobid_seff_Enformer_run=$(sbatch --dependency=afterany:$myjobid_Enformer_run --open-mode=append --output=$outfile_Enformer_run --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Enformer_run >> $outfile_Enformer_run")


#### interpreting_enformer_results #############################


type=$(echo "interpreting_enformer_results""_""$analysis")
outfile_interpreting_enformer_results=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_interpreting_enformer_results
echo -n "" > $outfile_interpreting_enformer_results
name_interpreting_enformer_results=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_interpreting_enformer_results=$(echo "$Rscripts_path""32_Interpret_Enformer_results.R")


Table_S6=$(echo "/group/soranzo/manuel.tardaguila/Motif_analysis/""Table_S6.tsv")
Enformer_result=$(echo "/group/soranzo/manuel.tardaguila/Enformer/Table_S6_variants/Enformer_results/variant_scores.all.tsv")
K562_features=$(echo "121_DNASE:K562,122_DNASE:K562,123_DNASE:K562,625_DNASE:K562")
HL60_features=$(echo "96_DNASE:HL-60")
context_initiator_TFs_subset=$(echo 'IRF3,IRF7,IRF1,IRF2,IRF9,IRF4,IRF8,RUNX2,RUNX1,RUNX3,FOSL1,FOS,FOSB,FOSL2,PU.1')
context_only_TFs_subset=$(echo 'FOXA1,FOXA2,FOXA3,FOXC,FOXP1,FOXP2,GATA3,MEF2C,LHX1,NKX6-1')
version_string=$(echo 'REF,ALT')
Intersect_SNP=$(echo 'YES')
TF_motifs_annotated=$(echo '/group/soranzo/manuel.tardaguila/Motif_analysis/Motifs_in_all_index_variants/TF_motifs_FINAL_K562_EXP.tsv')
Enformer_results=$(echo '/group/soranzo/manuel.tardaguila/Enformer/Table_S6_variants/""Enformer_summary_K562_HL60.tsv')



# myjobid_interpreting_enformer_results=$(sbatch --job-name=$name_interpreting_enformer_results --output=$outfile_interpreting_enformer_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_interpreting_enformer_results --Table_S6 $Table_S6 --Enformer_result $Enformer_result --K562_features $K562_features --HL60_features $HL60_features --context_initiator_TFs_subset $context_initiator_TFs_subset --context_only_TFs_subset $context_only_TFs_subset --version_string $version_string --Intersect_SNP $Intersect_SNP --TF_motifs_annotated $TF_motifs_annotated --Enformer_results $Enformer_results --type $type --out $output_dir")
# myjobid_seff_interpreting_enformer_results=$(sbatch --dependency=afterany:$myjobid_interpreting_enformer_results --open-mode=append --output=$outfile_interpreting_enformer_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_interpreting_enformer_results >> $outfile_interpreting_enformer_results")

myjobid_interpreting_enformer_results=$(sbatch --dependency=afterany:$myjobid_Enformer_run --job-name=$name_interpreting_enformer_results --output=$outfile_interpreting_enformer_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_interpreting_enformer_results --Table_S6 $Table_S6 --Enformer_result $Enformer_result --K562_features $K562_features --HL60_features $HL60_features --context_initiator_TFs_subset $context_initiator_TFs_subset --context_only_TFs_subset $context_only_TFs_subset --version_string $version_string --Intersect_SNP $Intersect_SNP --TF_motifs_annotated $TF_motifs_annotated --Enformer_results $Enformer_results --type $type --out $output_dir")
myjobid_seff_interpreting_enformer_results=$(sbatch --dependency=afterany:$myjobid_interpreting_enformer_results --open-mode=append --output=$outfile_interpreting_enformer_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_interpreting_enformer_results >> $outfile_interpreting_enformer_results")


#### enformer_graphs #############################


type=$(echo "enformer_graphs""_""$analysis")
outfile_enformer_graphs=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_enformer_graphs
echo -n "" > $outfile_enformer_graphs
name_enformer_graphs=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_enformer_graphs=$(echo "$Rscripts_path""33_Enformer_graphs.R")

Enformer_results=$(echo "/group/soranzo/manuel.tardaguila/Enformer/Table_S6_variants/""Enformer_summary_K562_HL60_TF_class.rds")
MPRA_results=$(echo "/group/soranzo/manuel.tardaguila/Enformer/Table_S5_MPRA_Results.tsv")


# myjobid_enformer_graphs=$(sbatch --job-name=$name_enformer_graphs --output=$outfile_enformer_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_enformer_graphs --Enformer_results $Enformer_results --MPRA_results $MPRA_results --type $type --out $output_dir")
# myjobid_seff_enformer_graphs=$(sbatch --dependency=afterany:$myjobid_enformer_graphs --open-mode=append --output=$outfile_enformer_graphs --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_enformer_graphs >> $outfile_enformer_graphs")

myjobid_enformer_graphs=$(sbatch --dependency=afterany:$myjobid_interpreting_enformer_results --job-name=$name_enformer_graphs --output=$outfile_enformer_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_enformer_graphs --Enformer_results $Enformer_results --MPRA_results $MPRA_results --type $type --out $output_dir")
myjobid_seff_enformer_graphs=$(sbatch --dependency=afterany:$myjobid_enformer_graphs --open-mode=append --output=$outfile_enformer_graphs --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_enformer_graphs >> $outfile_enformer_graphs")
