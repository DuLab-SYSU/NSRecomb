import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import CSS4_COLORS as mcolors
from collections import defaultdict


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
}


color_map = {
    'SARS-CoV-2': '#F80000',
    'RaTG13': '#4169E1',
    'PangolinGD': '#32CD32',
    'PangolinGX': '#EFEA3A',
    'BatSL': '#FF8C00',
    'RmYN02': '#9400D3',
    'SARS': '#FF69B4',
    'SARSr-CoVs': '#323232',
}

# color_map = {
#     'SARS-CoV-2': '#BC3C29',
#     'RmYN02': '#7876B1',
#     'RaTG13': '#0072B5',
#     'PangolinGD': '#20854E',
#     'PangolinGX': '#FFDC91',
#     'BatSL': '#E18727',
#     'SARS': '#EE4C97',
#     'SARSr-CoVs': '#323232',
# }


def SeqIdTitle(_in):
    data = {}
    with open(_in) as f:
        for line in f:
            k, v = line.strip().split('##')
            k2 = k.split('|')[1]
            data[k] = v
            data[k2] = v
    return data


def read_json(_in):
    with open(_in) as f:
        res = json.load(f)
    res_ = []
    for k, v in res.items():
        for v_i in v:
            v_i.append(int(k))
            res_.append(v_i)
    return res_


def get_ds_result(row):
    if (row[3] >= 0.8) and (row[4] > row[5]):
        return row[2]
    else:
        return row[1]


def get_bs_result(row):
    if row[3] >= 0.8:
        return row[2]
    else:
        return row[1]


def get_fragments(_in, df):
    fragments = []
    start = df.loc[0, 'site2']+251
    end = start
    status = df.loc[0, _in]

    if len(df) >= 10:
        for site, target in df[['site2', _in]].values:
            if target != status:
                fragments.append([status, start, end, end-start])
                status = target
                start = site+251
                end = site+251
            else:
                end = site+251
        else:
            fragments.append([status, start, end, end - start])
            fragments[0][1] = fragments[0][1] - 251
            fragments[-1][2] = fragments[-1][2] + 251
    else:
        fragments.append([status, start-251, start-251+226, 226])
    return fragments


def parse_result(_in, backbone, offset=251, length_threshold=199, log=False, plot=False):
    seq_id_title = SeqIdTitle('SeqIdTitle')
    res = read_json(_in)
    df = pd.DataFrame(
        res, columns=['query', 'backbone', 'subject', 'bs_value', 'ds1', 'ds2', 'site'])
    df['site2'] = df['site'] + offset
    df['result2'] = df.apply(get_ds_result, axis=1)
    df['subject_title'] = df['result2'].apply(
        lambda x: key_map.get(seq_id_title.get(x), 'SARSr-CoVs'))

    fragments = get_fragments('subject_title', df)
    fragments = [fragment for fragment in fragments if (
        fragment[0] != backbone) and (fragment[3] >= length_threshold)]

    if log:
        print(fragments)
    if plot:
        fig, host = plt.subplots()
        par1 = host.twinx()

        p1, = host.plot(df.site2, df.bs_value, label='Bootstrap support',
                        color=mcolors['grey'], alpha=0.8, linewidth=0.5, linestyle='solid')
        p2, = par1.plot(df.site2, df.ds1, label='dS(query-backbone)',
                        color=mcolors['grey'], alpha=0.8, linewidth=0.5, linestyle='dash')
        p3, = par1.plot(df.site2, df.ds2, label='dS(query-subject)',
                        color=mcolors['grey'], alpha=0.8, linewidth=0.5, linestyle='dotted')
        host.axhline(y=0.8, linestyle='--', color='r', alpha=0.6, linewidth=1)

        host.set_ylim([0, 1.5])
        host.set_ylabel('Bootstrap support')
        host.set_xlabel('Site')
        par1.set_ylabel('dS')

        ss = defaultdict(list)
        for id_, s_, e_, _ in fragments:
            ss[id_].append([s_, e_])

        lines = [p1, p2, p3]
        count = 0
        for k, v in ss.items():
            min_, max_ = zip(*v)
            tmp = host.hlines(y=[1.1]*len(min_), xmin=min_, xmax=max_,
                              label=k, color=color_map[k], linewidth=3)
            locals()['v_%s' % count] = tmp
            lines.append(locals()['v_%s' % count])
            count += 1

        host.legend(lines, [l_.get_label() for l_ in lines],
                    loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.show()
    return df, fragments


if __name__ == "__main__":
    backbone_list = ['EPI_ISL_412977', 'MN996532', 'MT121216', 'MT072864', 'MG772933']
    ss_ = ['SARS-CoV-2', 'RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL', 'SARSr-CoVs', 'SARS']
    job_list = ['SARS-CoV-2', 'RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX']
    seq_length = [[0, 21291, 3822, 228, 669],
                  [0, 21285, 3684, 228, 666],
                  [0, 21288, 3810, 231, 666],
                  [0, 21270, 3798, 228, 669],
                  [0, 21267, 3804, 228, 669]]

    fig, axs = plt.subplots(6, 1)
    axs[0].arrow(x=0, y=0, dx=21292, dy=0, head_width=0.1,
                 head_length=80, fc='k', ec='k')
    axs[0].arrow(x=21292, y=0, dx=3822, dy=0, head_width=0.1,
                 head_length=80, fc='k', ec='k')
    axs[0].arrow(x=25114, y=0, dx=228, dy=0, head_width=0.1,
                 head_length=80, fc='k', ec='k')
    axs[0].arrow(x=25339, y=0, dx=669, dy=0, head_width=0.1,
                 head_length=80, fc='k', ec='k')
    axs[0].arrow(x=26011, y=0, dx=1260, dy=0, head_width=0.1,
                 head_length=80, fc='k', ec='k')
    axs[0].set_ylim([-0.5, 0.5])
    axs[0].set_axis_off()

    for idx, job in enumerate(job_list):
        backbone = backbone_list[idx]
        cum_seq_length = np.cumsum(seq_length[idx])

        df1, fragments1 = parse_result("/home/zeng/python_work/nCov/result/%s.ORF1ab.result" % job,
                                       backbone, offset=cum_seq_length[0], length_threshold=6, plot=False, log=True)
        df2, fragments2 = parse_result("/home/zeng/python_work/nCov/result/%s.S.result" % job,
                                       backbone, offset=cum_seq_length[1], length_threshold=6, plot=False, log=True)
        df3, fragments3 = parse_result("/home/zeng/python_work/nCov/result/%s.E.result" % job,
                                       backbone, offset=cum_seq_length[2], length_threshold=6, plot=False, log=True)
        df4, fragments4 = parse_result("/home/zeng/python_work/nCov/result/%s.M.result" % job,
                                       backbone, offset=cum_seq_length[3], length_threshold=0, plot=False, log=True)
        df5, fragments5 = parse_result("/home/zeng/python_work/nCov/result/%s.N.result" % job,
                                       backbone, offset=cum_seq_length[4], length_threshold=6, plot=False, log=True)

        host = axs[idx+1]
        par1 = host.twinx()

        host.set_ylim([0, 1.5])
        if idx != 4:
            host.set_xticks([])
            host.set_xlabel('')
        else:
            host.set_xticks(range(0, 30001, 5000))
            host.set_xticklabels([0, 5000, 10000, 15000, 20000, 25000, 30000])
            host.set_xlabel('Site')
        if idx == 2:
            host.set_ylabel('Bootstrap support')
            par1.set_ylabel('dS')

        host.axhline(y=0.8, linestyle='--', color='r', alpha=0.8, linewidth=0.5)
        # host.scatter(x=cum_seq_length[1:], y=[1.1]*4, s=1)

        for df, fragments in [(df1, fragments1), (df2, fragments2), (df3, fragments3), (df4, fragments4), (df5, fragments5)]:
            p1, = host.plot(df.site2, df.bs_value, label='Bootstrap support',
                            color=mcolors['grey'], alpha=0.8, linewidth=0.5)
            p2, = par1.plot(df.site2, df.ds1, label='dS(query-backbone)',
                            color=mcolors['skyblue'], alpha=0.8, linewidth=0.5)
            p3, = par1.plot(df.site2, df.ds2, label='dS(query-subject)',
                            color=mcolors['bisque'], alpha=0.8, linewidth=0.5)

            ss = defaultdict(list)
            for id_, s_, e_, _ in fragments:
                ss[id_].append([s_, e_])

            count = 0
            for k, v in ss.items():
                min_, max_ = zip(*v)
                tmp = host.hlines(y=[1.1]*len(min_), xmin=min_, xmax=max_,
                                  label=k, color=color_map[k], linewidth=3)
                count += 1

        host.text(-1000, 1.55, 'Query: %s Backbone: %s' % (ss_[idx], ss_[idx+1]), )

    custom_lines = []
    for s_ in ss_:
        custom_lines.append(
            Line2D([0], [0], color=color_map[s_], lw=4, label=s_))
    # axs[1].legend(custom_lines, ss_, loc='upper left', bbox_to_anchor=(1.05, 1.2))
    plt.show()
