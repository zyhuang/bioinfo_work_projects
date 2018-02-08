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


def query(time_stamp, path, sub_level=0):

    path = path.rstrip('/')
    log_file = '_'.join(path.lstrip('/').split('/')[:4])
    if time_stamp != 'latest':
        log_file = 'out/{}.{}.log'.format(log_file, time_stamp)
    else:
        log_file_list = []
        for fname in os.listdir('out'):
            if log_file in fname:
                log_file_list.append(fname)
        log_file= 'out/' + sorted(log_file_list)[-1]

    dir_level = len(path.split('/')) + sub_level

    print(':: reading ' + log_file, file=sys.stderr)
    out_data = []
    for line in open(log_file):
        data = line.rstrip().split()
        if len(data) != 2:
            continue
        if path not in data[1]:
            continue
        level = len(data[1].split('/'))
        if level == dir_level:
            out_data.append(data)


    for d in sorted(out_data, key=lambda x: size_in_byte(x[0]), reverse=True):
        print(d[0], d[1], sep='\t')


if __name__ == '__main__':

    if len(sys.argv) != 4:
        print('python3  {}  [time_stamp|latest]  [path]  [sub_level]'
              .format(__file__))
        sys.exit(0)

    time_stamp, path, sub_level = sys.argv[1:4]
    sub_level = int(sub_level)
    query(time_stamp, path, sub_level)
