import json 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors 
import numpy as np 
from collections import Counter, defaultdict
from pprint import pprint
from scipy.stats import binom

key_map = {'RaTG13': 'RaTG13',  'pangolin/GD': 'pangolin/GD',  'bat_SL_CoVZXC21': 'bat_SL',  'bat_SL_CoVZC45': 'bat_SL',  'pangolin/GX': 'pangolin/GX',  'Longquan_140': 'Longquan_140',  'NA': 'NA',  'HKU3_1': 'HKU3', 'HKU3_2': 'HKU3', 'HKU3_4': 'HKU3', 'HKU3_6': 'HKU3', 'HKU3_3': 'HKU3', 'HKU3_7': 'HKU3', 'HKU3_8': 'HKU3', 'HKU3_13': 'HKU3', 'HKU3_5': 'HKU3', 'HKU3_9': 'HKU3', 'HKU3_11': 'HKU3', 'HKU3_12': 'HKU3', 'HKU3_10': 'HKU3'}  
labels = ['RaTG13', 'pangolin/GD', 'pangolin/GX', 'bat_SL', 'Longquan_140', 'HKU3', 'SARSr','NA'] 
color_map = {v: list(mcolors.TABLEAU_COLORS.values())[i] for i,v in enumerate(labels)} 


path = "/Users/jinfeng/Desktop/recombination-analysis/all_merged4.fas"



def load_MSA(_in):
    msa = defaultdict(dict)
    with open(_in) as f:
        for line in f:
            if line.startswith(">"):
                label = line[1:-1]
                msa[label] = ''
            else:
                msa[label] += line.strip().lower()
    return msa

msa = load_MSA(path)

# 'EPI_ISL_402125_SARS-CoV-2', 'MN996532_RaTG13', 'EPI_ISL_410721_pangolin/GD', 'EPI_ISL_410542_pangolin/GX', 'MG772933_bat_SL_CoVZC45', 'MG772934_bat_SL_CoVZXC21'
# sars_cov_2 = [(22621, 23016)]  
# ratg12 = [(1, 9108), (14440, 15216), (21145, 22476)] 
# pangolin_gd = [(12751, 16119), (17566, 20361), (22402, 24645), (25897, 27393)]

sars_cov_2 = msa['EPI_ISL_402125_SARS-CoV-2']
ratg13 = msa['MN996532_RaTG13']
pangolin_gd = msa['EPI_ISL_410721_pangolin/GD']
bat_sl = msa['MG772933_bat_SL_CoVZC45']

length = len(sars_cov_2)

def cal_p_value(_in, region): # _in (query, background, subject)
    q, b, s = _in
    start, end = region
    length = len(q)
    L, N, M = 0, 0, 0
    for i in range(length):
        if not q[i]==b[i]==s[i]:
            L += 1
            if start-1<=i<=end-1:
                N += 1
                if q[i] == s[i]:
                    M += 1
    p = len([i for i in range(length) if q[i] == s[i]]) / length
    p_value = 1*(L/N)*(1-binom.cdf(N-M, N, p))
    return L, N, M, p, p_value


L, N, M, p, p_value = cal_p_value((sars_cov_2, ratg13, pangolin_gd), (22621, 23016))
print(L, N, M, p, p_value)

L, N, M, p, p_value = cal_p_value((ratg13, pangolin_gd, bat_sl), (1, 9108))
print(L, N, M, p, p_value)

L, N, M, p, p_value = cal_p_value((ratg13, pangolin_gd, bat_sl), (21145, 22476))
print(L, N, M, p, p_value)

