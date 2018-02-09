'''Generate a list of unique pair of any items in a given list.

e.g. A,B,C,... => A,B A,C B,C ...
'''

import sys

def make_pair(list_name):
    '''generate a list of pair of items in an input list.

    Args:
        list_name (str):  input list name

    Return:
       (stdout): comma-separated pairs
    '''

    item_list = []
    for line in open(list_name):
        item_list.append(line.rstrip())

    item_list = sorted(item_list)
    for i,a in enumerate(item_list):
        for j in range(i, len(item_list)):
            b = item_list[j]
            print('{},{}'.format(a,b))


if __name__ == '__main__':

    list_name = sys.argv[1]
    make_pair(list_name)
