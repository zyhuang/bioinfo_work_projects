import os
import sys


def size_in_byte(size_str):

    power_hash = {
        'K': 1, 'M': 2, 'G': 3, 'T': 4,
    }
    if size_str[-1] not in power_hash:
        size = float(size_str)
    else:
        size = float(size_str[:-1]) * pow(1024, power_hash[size_str[-1]])
    return size


def stat_storage(timestamp):

    path_dict = {}
    for log_name in os.listdir('out'):
        if not log_name.endswith('.{}.log'.format(timestamp)):
            continue
        path = '/' + (log_name.rsplit('.', 2)[0].replace('_', '/'))
        for line in open('out/{}'.format(log_name)):
            data = line.rstrip().split()
            if data[1] == 'total':
                path_dict[path] = [data[0], size_in_byte(data[0])]

    print('', 'as of: {}'.format(timestamp), sep='\t')
    size2 = 0
    for path, size in sorted(path_dict.items(), key=lambda x: x[1][1], reverse=True):
        if 'snfs2' in path:
            print(path, size[0], size[1], sep='\t')
            size2 += size[1]
    print('snfs2 total:\t{0:.1f}T\n'.format(round(size2/pow(1024,4),2)))

    size4 = 0
    for path, size in sorted(path_dict.items(), key=lambda x: x[1][1], reverse=True):
        if 'snfs4' in path:
            print(path, size[0], size[1], sep='\t')
            size4 += size[1]
    print('snfs4 total:\t{0:.1f}T\n'.format(round(size4/pow(1024,4),2)))

    print('all total:\t{0:.1f}T'.format(round((size4+size2)/pow(1024,4),2)))

timestamp = sys.argv[1]
stat_storage(timestamp)
