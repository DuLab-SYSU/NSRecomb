import json
from collections import defaultdict
import matplotlib.pyplot as plt


def read_log(_in, bootstrap_threshold):
    with open(_in) as f:
        log = json.load(f)

    backbone_dict = defaultdict(list)
    hits_dict = defaultdict(list)
    bs_values = []
    for site, results in log.items():
        for result in results:
            query, backbone, hit, bs_value, ds1, ds2 = result
            # bootstrap value > 0.8 & dS(Q-B) > dS(Q-S)
            if bs_value > bootstrap_threshold and ds1 > ds2:
                hits_dict[hit].append(int(site))
            else:
                backbone_dict[backbone].append(int(site))
            bs_values.append((int(site), bs_value))
    bs_values_sorted = sorted(bs_values, key=lambda x: x[0])
    return backbone_dict, hits_dict, bs_values_sorted


def merge(sites, step):
    if len(sites) == 1:
        return [(sites[0], sites[0]+step-1)]
    else:
        segments = []
        start, end, last = sites[0], sites[0], sites[0]
        for site in sites[1:]:
            if site - last > step:
                segments.append((start, end+step-1))
                start = site
            else:
                end = site
            last = site
        else:
            # site = sites[-1]
            segments.append((start, last+step-1))
    return segments


def plot_result(bs_values_list, backbone_segments, hit_segments, outfig='tmp.jpg', interavtive_mode=False):
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    bs_sites, bs_values = zip(*bs_values_list)

    fig, ax = plt.subplots(1, 1, figsize=(10, 3))
    plt.rc('font', family='Arial', size=8)
    plt.subplots_adjust(bottom=0.2)

    ax.plot(bs_sites, bs_values, color='grey', alpha=0.5)
    ax.axhline(y=0.8, c='r', alpha=0.2, linestyle='--')

    count = 0
    for backbone, segments in backbone_segments.items():
        xmins, xmaxs = zip(*segments)
        ax.hlines(y=[1.5]*len(segments), xmin=xmins, xmax=xmaxs, label=backbone, colors=colors[count], linewidth=3)
        ax.text(x=max(bs_sites)+200, y=1.5, s=backbone, color=colors[count], verticalalignment='center')
    count += 1
    for hit, segments in hit_segments.items():
        xmins, xmaxs = zip(*segments)
        ax.hlines(y=[1.5-count*0.1]*len(segments), xmin=xmins, xmax=xmaxs, label=hit, colors=colors[count], linewidth=3)
        ax.text(x=max(bs_sites)+200, y=1.5-count*0.1, s=hit, color=colors[count], verticalalignment='center')
        count += 1

    ax.get_yaxis().set_visible(False)
    for spine in ["left", "top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.set_ybound(lower=-.1, upper=1.7)
    plt.xticks(fontproperties = 'Arial', size = 12)
    if interavtive_mode:
        plt.show()
    else:
        fig.savefig(outfig)


def parse_log(_in, step, gap_window, bootstrap_threshold, plot=False, interavtive_mode=False):
    backbone_dict, hits_dict, bs_values_list = read_log(_in, bootstrap_threshold)

    backbone_segments = dict()
    for backbone, results in backbone_dict.items():
        segments = merge(sorted(results), step)
        # print("Backbone: ", backbone, segments)
        backbone_segments[backbone] = segments
    hit_segments = dict()
    for hit, results in hits_dict.items():
        segments = merge(sorted(results), step)
        segments = [(start, end) for start, end in segments if end-start > step*gap_window]
        # print("Hits: ", hit, segments)
        if segments:
            hit_segments[hit] = segments

    if plot:
        outname = plot
        plot_result(bs_values_list, backbone_segments, hit_segments, outfig=outname, interavtive_mode=interavtive_mode)

    return hit_segments


if __name__ == '__main__':
    # parse_log("/home/zeng/python_work/nCov/result/tmp.log", step=15, plot=True, outfig='tmp.jpg')
    res = parse_log("/home/zeng/DulabWork/nCov/simulation/result/replicates_0_5/simu_replicates_0_5e-5_sequences00037_2.result", step=30, gap_window=0, bootstrap_threshold=0.8, plot='tmp', interavtive_mode=False)
    print(res)
