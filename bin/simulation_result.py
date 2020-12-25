import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_json(_in):
    with open(_in) as f:
        res = json.load(f)
    res_ = []
    for k, v in res.items():
        v.append(int(k))
        res_.append(v)
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
    start = 0
    end = 0
    status = df.loc[0, _in]
    for site, target in df[['site', _in]].values:
        if target != status:
            fragments.append([status, start+251, end+251, end-start])
            status = target
            start = site
            end = site
        else:
            end = site
    else:
        fragments.append([status, start+251, end+251, end - start])
    return fragments


def parse_result(_in, log=False, plot=False):
    res = read_json(_in)
    df = pd.DataFrame(res, columns=['query', 'backbone', 'subject', 'bs_value', 'ds1', 'ds2', 'site'])
    df['result1'] = df.apply(get_bs_result, axis=1)
    df['result2'] = df.apply(get_ds_result, axis=1)
    fragments1 = get_fragments('result1', df)
    fragments2 = get_fragments('result2', df)

    fragments1 = [fragments for fragments in fragments1 if (fragments[0] == 'S3') and (fragments[3] >= 102)]
    fragments2 = [fragments for fragments in fragments2 if (fragments[0] == 'S3') and (fragments[3] >= 102)]
    if log:
        print(fragments1)
        print(fragments2)
    if plot:
        df[['bs_value', 'ds1', 'ds2']].plot()
        plt.savefig('log___.pdf')
    return fragments1, fragments2


def main(j):
    result1 = []
    result2 = []
    for i in range(100):
        # print("*************** Replicates %s ******************" % i)
        r1, r2 = parse_result("../simulation/s%s/test%s_log_%s" % (j+1, j, i), log=False)

        total_length1 = 0
        for r_ in r1:
            total_length1 += r_[3]
        total_length2 = 0
        for r_ in r2:
            total_length2 += r_[3]
        result1.append(total_length1)
        result2.append(total_length2)
        # print("Total length: ", total_length1/1800, total_length2/1800)
    # print(result1)
    print(len([i for i in result1 if i != 0]) / 100, np.mean([i for i in result1 if i != 0]))

    # print(result2)
    print(len([i for i in result2 if i != 0]) / 100, np.mean([i for i in result2 if i != 0]))
    return result1, result2


if __name__ == "__main__":
    fragments = parse_result("../simulation/s3/test2_log_12", plot=True, log=True)
    # result1, result2 = main(1)
    # result1_ = [i for i in result1 if i != 0]
    # result2_ = [i for i in result2 if i != 0]

    # result3, result4 = main(2)
    # result3_ = [i for i in result3 if i != 0]
    # result4_ = [i for i in result4 if i != 0]

    # plt.figure()
    # plt.subplot(121)
    # plt.boxplot([result1_, result2_])
    # plt.grid(axis='y')
    # plt.xticks(ticks=[1, 2], labels=['Similarity-based', 'Synonymous distance\nadjusted'])
    # # plt.ylim([800, 2100])
    # plt.ylabel('Length of Recombination\nRegion', fontsize=12)
    # plt.title('Recombination', fontsize=12)

    # plt.table(cellText=['100%', '100%'], rowLabels=['Identification rate'], loc='bottom')

    # plt.subplot(122)
    # plt.boxplot([result3_, result4_])
    # plt.grid(axis='y')
    # plt.xticks(ticks=[1, 2], labels=['Similarity-based', 'Synonymous distance\nadjusted'])
    # plt.ylim([0, 1000])
    # plt.title('Localized Convergent Evolution', fontsize=12)
    # plt.margins(y=0.05)
    # plt.table(cellText=['70%', '11%'], rowLabels=['Identification rate'], loc='bottom')
    
    # plt.show()
