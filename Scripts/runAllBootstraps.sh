#BSUB -n 2
#BSUB -q general
#BSUB -W 4:00

#BSUB -o Joboutputs/parallelBootstrap.%I.%J.out
#BSUB -e Joboutputs/parallelBootstrap.%I.%J.err
#BSUB -J bootstrap[1-100]

echo "Job: ${LSB_JOBINDEX}"
source activate bio3
python parallelBootstrap.py /research/gutsybugs/HLA/Data/pops_A-B-DQB1_min95.pcl maj_ABDQB1_min95 ${LSB_JOBINDEX}