import re
import json
from os.path import abspath, dirname
import sys


def load_json(fname):
    return json.loads(open(fname).read())


field = None
json_name = sys.argv[1]

if len(sys.argv) == 3:
    field = sys.argv[2]

data = load_json(json_name)

if not field:

    print(json.dumps(data, sort_keys=True, indent=4))

else:
    value = data
    for key in field.split('.'):
        index = -1
        if '[' in key and ']' in key:
            key, index = key.replace(']','').split('[')
            index = int(index)
        if index >= 0:
            value = value[key][index]
        else:
            value = value[key]

    print(json.dumps(value, sort_keys=True, indent=4))
