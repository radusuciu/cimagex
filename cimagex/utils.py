from itertools import tee, filterfalse

def partition(pred, iterable):
    'Use a predicate to partition entries into true entries and false entries'
    # partition(is_odd, range(10)) --> 1 3 5 7 9  and 0 2 4 6 8
    t1, t2 = tee(iterable)
    return list(filter(pred, t2)), list(filterfalse(pred, t1))
