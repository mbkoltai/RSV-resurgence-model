# master start file
#!/usr/bin/bash 
scp batch_run_files/start_batches.sh .
scp batch_run_files/batch*.sh .
module load R/3.6.3
sh start_batches.sh
rm batch*.sh
rm start_batches.sh
