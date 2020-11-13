#!/usr/bin/env python3

import subprocess
import os

RESULT_TAG = "[[RESULT]]"

def model_call(sampleData):
    param = sampleData["Parameters"]
    beta = param[0]
    sigma = param[1]

    p = subprocess.Popen(["./model_run.py", str(beta), str(sigma)],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    lines = p.stdout.readlines()
    lines = [l.decode("utf-8") for l in lines]

    # write the output to log
    with open('model.out', 'a') as f:
        f.write(''.join(lines))

    lines_err = p.stderr.readlines()
    lines_err = [l.decode("utf-8") for l in lines_err]
    with open('model.err', 'a') as f:
        f.write(''.join(lines_err))

    ires = None
    for i, line in enumerate(lines):
        if line.strip() == RESULT_TAG:
            ires = i
    assert ires, "{} not found in process output.\nThe output was:{:}\n".format(
        RESULT_TAG, ''.join(lines))
    d = dict()
    exec('res = {' + ''.join(lines[ires + 1:]) + '}', d)
    res = d['res']
    for k in res:
        sampleData[k] = res[k]


def example():
    s = {"Parameters": [1, 2]}
    model_call(s)
    print(s)

if __name__ == "__main__":
    example()
