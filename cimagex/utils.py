from itertools import tee, filterfalse
import numpy as np

np.seterr(all='raise')

def partition(pred, iterable):
    'Use a predicate to partition entries into true entries and false entries'
    # partition(is_odd, range(10)) --> 1 3 5 7 9  and 0 2 4 6 8
    t1, t2 = tee(iterable)
    return list(filter(pred, t2)), list(filterfalse(pred, t1))


def get_modified_z_scores(values):
    non_zero_values = [y for y in values if y != 0]
    median = np.median(non_zero_values)
    mad = np.median([np.abs(y - median) for y in non_zero_values])
    modified_z_scores = [(0.6745 * (y - median) / mad) if (y != 0 and mad != 0) else 0 for y in values]
    return modified_z_scores

def get_iqr_bounds(values):
    values = sorted([v for v in values if v != 0])
    q1, q3 = np.percentile(values, [25, 75])
    lower = q1 - (1.5 * q1)
    upper = q3 + (1.5 * q3)
    return (lower, upper)
