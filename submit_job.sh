#BSUB -L /bin/bash
#BSUB -W 1:00
#BSUB -n 1
#BSUB -e ~/hpc_log/%J.err
#BSUB -o ~/hpc_log/%J.out

module load gcc/4.8.2 cdo/1.6.4
module load new mesa/12.0.6 rstudio/0.99.669 netcdf/4.3.2

echo echo test | ./rungpmodel