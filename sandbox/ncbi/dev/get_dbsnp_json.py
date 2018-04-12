'''This program downloads dbsnp json data using NCBI API.
'''

import sys
import json
import requests


rsid = sys.argv[1]

if rsid.startswith('rs'):
    rsid = rsid[2:]

url = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{}'.format(rsid)
out_json = '{}.json'.format(rsid)

print('get ' + url)
r = requests.get(url)
data = json.loads(r.text)

print('write ' + out_json)
fout = open(out_json, 'w')
print(json.dumps(data, sort_keys=True, indent=4), file=fout)
fout.close()
