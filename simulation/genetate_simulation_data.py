import os
import sys
from pathlib import Path
import subprocess


RECODON_PATH = Path.home().joinpath('opt', 'Recodon-1.6.0', 'exe')
# Population parameters
RECOMBINATION_RATES = ['0', '0.25e-6', '1e-6', '4e-6', '16e-6']
MUTATION_RATES = ['0.25e-5', '1.25e-5', '2.5e-5', '5e-5']
# Evolution model parameters
PI_F3x4 = '12 0.3 0.2 0.4 0.1 0.2 0.2 0.3 0.3 0.2 0.2 0.1 0.5'
KAPPA = '2.0'
OMEGA = '1 0.63'

# r = RECOMBINATION_RATES[1]
# mua = MUTATION_RATES[0]


rootPath = Path.cwd()

for i, r in enumerate(RECOMBINATION_RATES):
    for j, mua in enumerate(MUTATION_RATES):
        child_dir = "replicates_%s_%s" % (r, mua)
        if r == '0':
            current_dir = rootPath.joinpath('norecombinate', child_dir)
        else:
            current_dir = rootPath.joinpath('recombinate', child_dir)
        if not current_dir.exists():
            current_dir.mkdir(parents=True)

        print('recodon1.6.0_Linux -n100 -s10 -l1000 -e1000 -_1 -p0 -r%s -q1 10 0.00 -u%s -f%s -t%s -m%s -bseq%s%s -z -*2 -jtree%s%s -dbreakpoint%s%s -#123456'
            % (r, mua, PI_F3x4, KAPPA, OMEGA, i, j, i, j, i, j))

        result = subprocess.run(
            'recodon1.6.0_Linux -n100 -s10 -l1000 -e1000 -_1 -p0 -r%s -q1 10 0.00 -u%s -f%s -t%s -m%s -bseq%s%s -z -*2 -jtree%s%s -dbreakpoint%s%s -#123456'
            % (r, mua, PI_F3x4, KAPPA, OMEGA, i, j, i, j, i, j),
            shell=True,
            # capture_output=True,
            text=True,
            env={'PATH': RECODON_PATH},
            cwd=current_dir
        )
        
        # os.chdir(rootPath)
        # print(result.stdout)
        # print(result.stderr)
