import json
import os
import matplotlib.pyplot as plt
from collections import defaultdict

key_map = {
    'RaTG13': 'RaTG13',
    'pangolin/GD': 'PangolinGD',
    'bat_SL_CoVZXC21': 'BatSL',
    'bat_SL_CoVZC45': 'BatSL',
    'pangolin/GX': 'PangolinGX',
    'NA': 'NA'
}
color_list = [
    '#EB79CB', '#323232', '#F80000', '#4169E1', '#32CD32', '#FFFF00', '#FFA500'
]
lineages = [
    'SARS', 'SARSr-CoVs', 'SARS-CoV-2', 'RaTG13', 'PangolinGD', 'PangolinGX',
    'BatSL'
]

color_map = {a: b for a, b in zip(lineages, color_list)}

dirpath = '/Users/jinfeng/Desktop/evolution_and_diversity/blast_analysis'

tasks = [
    '2019-ncov-cds.0.8.json', 'ratg13-cds.0.8.json',
    'pangolin-gd-cds.0.8.json', 'pangolin-gx-cds.0.8.json'
]
querys = [
    'SARS-CoV-2', 'RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL',
    'SARSr-CoVs'
]
backbones = ['RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL']


def read_data(_in):
    # read data
    with open(_in) as f:
        contents = json.loads(f.read())
    results = {
        int(k): set(map(lambda x: key_map.get(x[1], 'SARSr-CoVs'), v))
        for k, v in contents.items()
    }
    x_ = [int(k) for k in contents.keys()]
    y_ = [sum([v_i[2] for v_i in v]) / len(v) for v in contents.values()]
    y_ = list(map(lambda x: 0.65 if x == 0 else x, y_))
    # unify by source virus
    source_by_virus = defaultdict(list)
    for k, v in results.items():
        for v_i in v:
            source_by_virus[v_i].append(k)
    return source_by_virus, x_, y_


def fragment(hits, inter_length=3, min_length=3):
    fragments = []
    start = hits[0]
    for i in range(len(hits) - 1):
        if hits[i + 1] > hits[i] + inter_length:
            end = hits[i] + inter_length - 1
            fragments.append((start, end, end - start + 1))
            start = hits[i + 1]
    else:
        end = hits[-1] + inter_length - 1
        fragments.append((start, end, end - start + 1))
    return [fragment for fragment in fragments if fragment[2] > min_length]


def plot_hlines(ax, task, query, **kwargs):
    order = tasks.index(task)
    source_by_virus, x_, y_ = read_data(os.path.join(dirpath, task))
    fragments_by_virus = {
        k: fragment(v, **kwargs)
        for k, v in source_by_virus.items()
    }
    # bootstrap plot
    ax.plot(x_, y_, alpha=0.8, lw=1, color='#808080', antialiased=True)
    # threshod line 0.8
    ax.axhline(y=0.8, alpha=0.8, lw=1, ls='--', c='r')

    yvalue = 1.05
    for k in querys:
        v = fragments_by_virus.get(k)
        if v:
            # print(k, v)
            xmins, xmaxs, _ = zip(*v)
            ax.hlines(y=[yvalue] * len(v),
                      xmin=xmins,
                      xmax=xmaxs,
                      colors=color_map[k],
                      lw=4,
                      label=k)
            yvalue += 0.05
            # ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., title='Query: %s' % query, title_fontsize=12)
    ax.set_xlim([-500, 28000])
    ax.set_ylim([0, 1.3])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=8)
    ax.text(0.1, 0.1, "Query: %s, Backbone: %s" % (querys[order], backbones[order]), fontsize=10, fontfamily='Arial')
    # ax.set_ylabel('bootstrap value', fontsize=8)
    # secax = ax.secondary_yaxis('right')
    # secax.set_yticks([])
    # secax.set_ylabel('%s' % query, fontsize=1)


def main(**kwargs):
    # tasks = ['sars-cov-2_80_500_cds.json', 'ratg13_80_500_cds.json', 'pangolin-gd_80_500_cds.json', 'pangolin-gx_80_500_cds.json']

    n = len(tasks)

    fig, axes = plt.subplots(n + 1, 1, figsize=(10, 6), sharex=True, dpi=300)
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=1, hspace=0)

    # 1-21360 21361-25224 25225-25455 25456-26124 26125-27393
    ax = axes[0]
    ax.set_xlim([0, 28000])
    ax.set_ylim([0, 0.4])
    ax.set_axis_off()
    ax.text((1+13203 - 250 - 80)/2, 0.11, 'ORF1a', fontsize=8, fontfamily='Arial')
    ax.arrow(1,
             0.1,
             13203 - 250 - 80,
             0,
             head_width=0.02,
             head_length=80,
             fc='k',
             ec='k')
    ax.text((13200 - 250 + 21292 - 250)/2, 0.11, 'ORF1b', fontsize=8, fontfamily='Arial')
    ax.arrow(13200 - 250,
             0.1,
             8087 - 80,
             0,
             head_width=0.02,
             head_length=80,
             fc='k',
             ec='k')
    ax.text((21292 - 250 + 25114 - 250)/2, 0.06, 'S', fontsize=8, fontfamily='Arial')
    ax.arrow(21292 - 250,
             0.05,
             3822 - 80,
             0,
             head_width=0.02,
             head_length=80,
             fc='k',
             ec='k')
    ax.text(25114 - 250, 0.11, 'E', fontsize=8, fontfamily='Arial')
    ax.arrow(25114 - 250,
             0.1,
             228 - 80,
             0,
             head_width=0.02,
             head_length=80,
             fc='k',
             ec='k')
    ax.text(25432 - 250, 0.06, 'M', fontsize=8, fontfamily='Arial')
    ax.arrow(25342 - 250,
             0.05,
             669 - 80,
             0,
             head_width=0.02,
             head_length=80,
             fc='k',
             ec='k')
    ax.text(26011 - 250, 0.11, 'N', fontsize=8, fontfamily='Arial')
    ax.arrow(26011 - 250,
             0.1,
             1260 - 250 - 80,
             0,
             head_width=0.02,
             head_length=80,
             fc='k',
             ec='k')

    # axes[0].set_title("Window_size=500, Step=3, Bootstrap replicates=500, Bootstrap cutoff=0.8")
    for i, task in enumerate(tasks):
        plot_hlines(axes[i + 1], task, querys[i], **kwargs)
    # plt.show()
    plt.xlabel("Genomic position", fontsize=10, fontfamily='Arial')
    # plt.ylabel("Bootstrp value", fontsize=8)
    plt.xticks(fontsize=8, fontfamily='Arial')
    fig.savefig('/Users/jinfeng/Desktop/evolution_and_diversity/figure/figure1b.svg')


main(inter_length=3, min_length=100)
