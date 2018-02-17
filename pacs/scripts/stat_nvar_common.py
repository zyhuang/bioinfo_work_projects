import sys
import json


def stat_nvar_common(binkey):

    counter = {}
    for line in open('af_prov_list/{}.af.list'.format(binkey)):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)
        for prov in vardata:
            if prov not in counter:
                counter[prov] = 0
            if vardata[prov]['ns_all'] > 0:
                counter[prov] += 1

    print(json.dumps(counter, sort_keys=True))


def stat_nvar_common_stdin():


    counter = {}
    for line in sys.stdin:
        data = json.loads(line.rstrip())
        for prov in data:
            if prov not in counter:
                counter[prov] = 0
            counter[prov] += data[prov]

    for prov, nvar in sorted(counter.items(), key=lambda x: x[1], reverse=True):
        print(prov, nvar, sep='\t')



if len(sys.argv) == 2:

    binkey = sys.argv[1]
    stat_nvar_common(binkey)

else:
    # cat log/stat_nvar_common.o1779460-* | grep -v ^% | python3 stat_nvar_common2.py
    stat_nvar_common_stdin()
