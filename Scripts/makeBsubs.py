import subprocess
import glob, os
bsub = '''#BSUB -n 2
#BSUB -q general
#BSUB -W 4:00

#BSUB -o Joboutputs/%s.%%I.%%J.out
#BSUB -e Joboutputs/%s.%%I.%%J.err
#BSUB -J bootstrap[1-100]
source activate bio3
'''
cmd = 'python parallelBootstrap.py %s %s ${LSB_JOBINDEX}\n'
for pcl in glob.glob('/research/btc_bioinformatic/operations/HLA/FreqRT/Data/Pops/pops*.pcl'):
    basename = os.path.splitext(os.path.basename(pcl))[0]
    with open(f'bsub_{basename}.sh' , 'rb') as w:
        w.write(bsub % (basename, basename))
        w.write(cmd % (pcl, basename))
        