import json
import os
import numpy as np
import matplotlib.pyplot as plt


def read_json(_in):
    with open(_in) as f:
        res = json.load(f)
    res_ = []
    for k, v in res.items():
        v.append(int(k))
        res_.append(v)
    return sorted(res_, key=lambda x: x[6])


if __name__ == "__main__":
    print(os.getcwd())
    res = read_json("../simulation/s2/test1_log_0")
    res_ = np.array(res)
    plt.plot(res_[:, -1], res_[:, -4])
    plt.show()