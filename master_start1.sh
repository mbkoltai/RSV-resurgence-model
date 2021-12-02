# master start file 
#!/usr/bin/bash 
Rscript fcns/write_run_file.R 64 1084 25 4 repo_data/partable_filtered_AR_seasconc.csv NOSAVE sep_qsub_files 2018-09-01 8
scp batch_run_files/start_batches.sh . 
scp batch_run_files/batch*.sh . 
module load R/3.6.3
sh start_batches.sh 
rm batch*.sh 
rm start_batches.sh
