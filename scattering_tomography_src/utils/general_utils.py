import operator


def Extract(lst, ind=0):
    return list(map(operator.itemgetter(ind), lst))
