import sys
from pathlib import Path
from collections import defaultdict

# add package to python module path
p = Path.cwd()
root_dir = p.parent
bin_dir = root_dir.joinpath("source")
sys.path.append(str(bin_dir))

import Merge


def read_log(log_file, step, gap_window, bootstrap_threshold):
    log_file = Path(log_file)
    _, _, _, _, replicate_i, sample_i  = log_file.stem.split("_")
    replicate_id = int(replicate_i.lstrip('sequences'))
    sample_id = int(sample_i)
    hit_segments = Merge.parse_log(log_file, step, gap_window, bootstrap_threshold, plot=False)
    return replicate_id, sample_id, hit_segments


def judge_recomb(bs_value, n_gap, result_dir):
    result_dir = Path(result_dir)
    s1_dirs = [x for x in result_dir.iterdir() if x.stem.startswith("replicate")]
    s2_dirs = [x for x in result_dir.iterdir() if x.stem.startswith("convergent")]

    res = []
    for dirs in [s1_dirs, s2_dirs]:
        print("Task")
        res_i = []
        for dir_i in dirs:
            print("="*100 + ">")
            title, r, mu = dir_i.name.split("_")
            print("recombination rate: %s, mutataion rate: %s" % (r, mu))
            logfiles = dir_i.glob("simu_*")
            res_i_j = defaultdict(list)
            for logfile in logfiles:
                replicate_id, sample_id, hit_segments = read_log(logfile, 30, n_gap, bs_value)
                if hit_segments:
                    res_i_j[replicate_id].append(1)
                    print(logfile.name)
                else:
                    res_i_j[replicate_id].append(0)
            count = 0
            for k, v in res_i_j.items():
                if sum(v) != 0:
                    count += 1
            res_i.append((r, mu, count))
            print(count)
            print("<" + "="*100, "\n")
        res.append(res_i)

    return res



result_dir = Path.cwd().joinpath('result')
res = judge_recomb(0.95, 7, str(result_dir))
for res_i in res:
    print("Task")
    print('%17s%15s%13s' % ('Reombination rate', 'Mutation rate', 'Probability'))
    res_sorted = sorted(res_i)
    for r, u, count in res_sorted:
        print("%17s%15s%13d" % (r, u, count))
