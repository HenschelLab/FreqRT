"""This script generates small bsub (IBM's LSF) job queue sub mission scripts, and launches them using 
subprocess.call. The scipts call the python program parallelBootstrap with three parameters:
(it uses job arrays)
    1. the pcl file containing all the population set data
    2. the basename
    3. the job array index
    
    Currently hardcoded for 100 bootstraps, change simply to 1-1000 in #BSUB -J directive"""

import subprocess
import glob, os
import time
bsub = '''#BSUB -n 2
#BSUB -q general
#BSUB -W 4:00

#BSUB -o Joboutputs/%s.%%I.%%J.out
#BSUB -e Joboutputs/%s.%%I.%%J.err
#BSUB -J %s[1-100]
source activate bio3
'''
cmd = 'python parallelBootstrap.py %s %s ${LSB_JOBINDEX}\n'
for pcl in glob.glob('/research/btc_bioinformatic/operations/HLA/FreqRT/Data/Pops/pops*.pcl')[1:]: ## first was done manually
    basename = os.path.splitext(os.path.basename(pcl))[0]
    print(basename)
    with open(f'bsub_{basename}.sh' , 'w') as w:
        w.write(bsub % (basename, basename, basename))
        w.write(cmd % (pcl, basename))
        #subprocess.call(f'bsub < bsub_{basename}.sh', shell=True)
        time.sleep(0.5)
    print(f'bsub < bsub_{basename}.sh')
        

