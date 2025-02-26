from statistics import mean
import numpy as np


def meanSize(x_axis_distribution, numFrac):
    d4 = 0
    d3 = 0
    for (diam, ni) in zip(x_axis_distribution, numFrac):
        d3 += ni * (diam**3)
        d4 += ni * (diam**4)

    meanSize = d4/d3
    return meanSize

def maxSize(x_axis_distribution, y_axis_pdf):
    mask = np.ones(np.size(y_axis_pdf), dtype=bool)
    for i in range(np.size(y_axis_pdf)):
        if y_axis_pdf[i] == 0:
            mask[i] = False
    x_axis_distribution_clean = x_axis_distribution[mask,...]
    maxSize = x_axis_distribution_clean[-1]
    return maxSize
