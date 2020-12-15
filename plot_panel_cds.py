import json
from collections import defaultdict
import matplotlib.pyplot as plt
import pprint


key_map = {
    'RmYN02': 'RmYN02',
    'RaTG13': 'RaTG13',
    'Pangolin/GD/MP789': 'PangolinGD',
    'PangolinGD': 'PangolinGD',
    'bat_SL_CoVZXC21': 'BatSL',
    'bat_SL_CoVZC45': 'BatSL',
    'BatSL': 'BatSL',
    'Pangolin/GX/P1E': 'PangolinGX',
    'Pangolin/GX/P5L': 'PangolinGX',
    'Pangolin/GX/P2V': 'PangolinGX',
    'Pangolin/GX/P4L': 'PangolinGX',
    'Pangolin/GX/P5E': 'PangolinGX',
    'PangolinGX': 'PangolinGX',
    'NA': 'NA'
}

# color_map = {
#     'SARS-CoV-2': '#F80000', 
#     'RaTG13' :'#4169E1', 
#     'PangolinGD': '#32CD32', 
#     'PangolinGX': '#FFFF00',
#     'BatSL': '#FF8C00', 
#     'RmYN02': '#9400D3',
#     'SARS': '#FF69B4', 
#     'SARSr-CoVs': '#323232', 
# }

color_map = {
    'SARS-CoV-2': '#BC3C29', 
    'RaTG13' :'#0072B5', 
    'PangolinGD': '#20854E', 
    'PangolinGX': '#FFDC91',
    'BatSL': '#E18727', 
    'RmYN02': '#7876B1',
    'SARS': '#EE4C97CC', 
    'SARSr-CoVs': '#323232', 
}


subjects_order = ['SARS-CoV-2', 'RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL', 'SARSr-CoVs']


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


def read_data(_in, **kwargs):
    # read data
    with open(_in) as f:
        contents = json.loads(f.read())
    
    results = {
        int(k): list(set(map(lambda x: key_map.get(x[1], 'SARSr-CoVs'), v)))
        for k, v in contents.items()
    }
    x_y = [(int(k), sum([float(v_i[2]) for v_i in v]) / len(v)) for k, v in contents.items()]
    # unify by source virus
    source_by_virus = defaultdict(list)
    for k, v in results.items():
        for v_i in v:
            source_by_virus[v_i].append(k)

    fragments_by_virus = {
        k: fragment(v, **kwargs)
        for k, v in source_by_virus.items()
    }
    return fragments_by_virus, x_y


def merge_cds(task):
    cds = ['ORF1ab', 'S', 'M', 'N']
    fragments_by_virus_merged = defaultdict(list)
    x_y_merged = []
    offset = 0
    for cds_i in cds:
        fragments_by_virus, x_y  = read_data("/home/zeng/python_work/bioinfo/nCov/recombination_task/task_yn02/%s/%s_%s.501.json" % (task, task,cds_i), inter_length=3, min_length=10)
        print("CDS:", cds_i)
        pprint.pprint(fragments_by_virus)
        x_y = sorted(x_y, key=lambda x:x[0])
        length = x_y[-1][0]
        x_y_merged.extend([(x_+offset, y_) for x_, y_ in x_y])
        f = lambda x: (x[0]+offset, x[1]+offset, x[2])
        for k, v in fragments_by_virus.items():
            fragments_by_virus_merged[k].extend(list(map(f, v)))

        offset += length
    return fragments_by_virus_merged, x_y_merged


def cds(ax, start, length, y, y2, label):
    ax.set_axis_off()
    ax.text(start+length*0.5-200, y2, label, fontsize=8, fontfamily='Arial')
    ax.arrow(start, y, length-80, 0, head_width=0.02, head_length=80, fc='k', ec='k')
        

def plot_panel(tasks=['SARS-CoV-2', 'RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX'], 
               backbones = ['RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL']):
    
    n_row = len(tasks)
    fig, axes = plt.subplots(n_row+1, 1, figsize=(5.5, 4), sharex=True, dpi=300)
    

    cds(axes[0], 1, 21291-500, 0.3, 0.33, 'ORF1ab')
    cds(axes[0], 21291-500, 3822-500, 0.2, 0.01, 'S')
    cds(axes[0], 21291+3822-1000, 228-100, 0.3, 0.33, 'E')
    cds(axes[0], 21291+3822+228-1000, 441-250, 0.2, 0.01, 'M')
    cds(axes[0], 21291+3822+228+441-1350, 1260-500, 0.3, 0.33, 'N')
    axes[0].set_ylim(0, 1)
    for i,task in enumerate(tasks):
        print(task)
        fragments_by_virus, x_y = merge_cds(task)
        x_, y_ = zip(*x_y)
        # bootstrap plot
        axes[i+1].plot(x_, y_, alpha=0.8, lw=0.3, color='#808080', antialiased=True)
        # threshod line 0.8
        axes[i+1].axhline(y=0.8, alpha=0.8, lw=0.5, ls='--', c='r')

        yvalue = 1.05
        for k in subjects_order:
            v = fragments_by_virus.get(k)
            if v:
                xmins, xmaxs, _ = zip(*v)
                axes[i+1].hlines(y=[yvalue] * len(v),
                        xmin=xmins,
                        xmax=xmaxs,
                        colors=color_map[k],
                        lw=3,
                        label=k)
                yvalue += 0.05
                # axes[i].legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        axes[i+1].set_xlim([-500, 26000])
        axes[i+1].set_ylim([0, 1.3])
        axes[i+1].set_yticks([0, 0.5, 1])
        axes[i+1].set_yticklabels([0, 0.5, 1], fontsize=8)
        axes[i+1].text(0.1, 0.1, "Query: %s, Backbone: %s" % (task, backbones[i]), fontsize=8, font='Arial')

    plt.xlabel("Genomic position", fontsize=10, font='Arial')
    #plt.ylabel("Bootstrp value", fontsize=8)
    plt.xticks(fontsize=8)

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=1, hspace=0)
    # fig.savefig("/home/zeng/Desktop/recombination_nCoV/recombination_task/Fig1_b.jpg")
    plt.show()


plot_panel()
