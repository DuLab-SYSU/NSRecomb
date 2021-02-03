from pathlib import Path
import subprocess
import re
from tempfile import NamedTemporaryFile


# the path for 3seq P value table
PTABLE = "/home/zeng/opt/3seq/PVT.3SEQ.2017.700"
PY3SEQ_PATH = '/home/zeng/opt/3seq'


def py3seq(_in):
    id_ = NamedTemporaryFile(mode='w+t', delete=True)

    client = subprocess.run(
        "3seq -full %s -d -ptable '%s' -id %s" % (_in, PTABLE, id_.name),
        shell=True,
        input=b"Y",
        capture_output=True,
        env={"PATH": PY3SEQ_PATH}
    )

    out = client.stdout.decode('utf-8')
    err = client.stderr.decode('utf-8')
    id_.close()
    try:
        res = int(out.split('\n')[37].split('\t')[1])
    except:
        print(_in)
        res=0
    return res


results = []

tasks = ['convergent', 'norecombinate', 'recombinate']
simu_dir = Path("/home/zeng/DulabWork/nCov/simulation/")

simu_dir = simu_dir / tasks[2]
simu_dirs = [x for x in simu_dir.iterdir() if x.is_dir()]

for simu_dir_i in simu_dirs:
    _, r, mua = simu_dir_i.name.split("_")
    replicates = simu_dir_i.joinpath('Results').glob("seque*")
    count = 0
    for replicate in replicates:
        if py3seq(replicate):
            count+=1
    results.append((r, mua, count))

# print(sorted(results))

print("%17s%15s%13s" % ('Recombination Rate', 'Mutataion Rate', 'Probability'))
for r, mua, count in sorted(results):
    print("%17s%15s%13s" % (r, mua, count))