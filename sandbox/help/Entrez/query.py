from Bio import Entrez
import sys
import json



def query(term_id, email="david.henke@bcm.edu"):

    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=term_id, retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    fout = open(term_id + '.json', 'w')
    print(json.dumps(record, sort_keys=True, indent=2), file=fout)
    fout.close()

    return record


if __name__ == '__main__':

    term_id = sys.argv[1]
    query(term_id)
