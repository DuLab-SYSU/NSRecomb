import json
import matplotlib as plt


def read_json(_in):
    with open(_in) as f:
        res = json.load(f)
    return res


def unify(data):
    idx = []
    res = []
    for k, v in data.items():
        for v_i in v:
            idx.append(k)
            res.append(v)

    querys, backbones, subjects, bs_values, dS1, dS2 = *zip(res)
    return idx, res


if __name__ == "__main__":
    data = read_json('log')
    unify(data)
