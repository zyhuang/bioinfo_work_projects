'''This progrom parses the pubmed xml data and save the results in a line-based
key-value json list file.

The XML data is downloaded from the following FTP
ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/
ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/

The extracted information about each article are

article:
- pmid
- doi
- title
- language
- date
    - day
    - month
    - year
- journal
    - title
    - volume
    - issue
    - pages
    - country
    - language
- abstract
- author_list (list of dict)
    - firstname
    - lastname
    - affiliations (list of str)
- keyword_list (list of str)
- ref_list (list of str)

'''

import sys
import json
from lxml import etree
import gzip


def get_value(elem, xpath_str, default=None):

    try:
        value = elem.xpath(xpath_str)[0].text
    except:
        value = default
    return value


def get_authors(elem):

    author_list = []
    for el in elem.xpath('.//AuthorList/Author'):
        lastname = get_value(el, './LastName')
        firstname = get_value(el, './ForeName')
        affiliations = []
        for af in el.xpath('./AffiliationInfo/Affiliation'):
            affiliations.append(af.text)
        author_list.append(
            {'firstname': firstname, 'lastname': lastname,
             'affiliations': affiliations}
        )
    return author_list


def get_keywords(elem):

    keyword_list = []
    for el in elem.xpath('.//KeywordList/Keyword'):
        keyword_list.append(el.text)
    return keyword_list


def get_references(elem):

    ref_list = []
    for el in elem.xpath('.//CommentsCorrections'):
        ref_list.append(get_value(el, './PMID'))
    return ref_list


def get_abstract(elem):

    try:
        abstract = ''.join(elem.xpath('.//AbstractText')[0].itertext()).strip()
    except:
        abstract = None
    return abstract


def get_article_info(article_el):

    data = dict()
    extract = {
        'title': './/ArticleTitle',
        'pmid': './/PMID',
        'doi': './/ELocationID',
    }

    extract_date = {
        'day': './/JournalIssue/PubDate/Day',
        'month': './/JournalIssue/PubDate/Month',
        'year': './/JournalIssue/PubDate/Year',
    }

    extract_journal = {
        'title': './/Journal/Title',
        'volume': './/JournalIssue/Volume',
        'issue': './/JournalIssue/Issue',
        'pages': '..//Pagination/MedlinePgn',
        'country': './/MedlineJournalInfo/Country',
        'language': './/Language',
    }

    for k, v in extract.items():
        data[k] = get_value(article_el, v)

    data['date'] = {}
    for k, v in extract_date.items():
        data['date'][k] = get_value(article_el, v)

    data['journal'] = {}
    for k, v in extract_journal.items():
        data['journal'][k] = get_value(article_el, v)

    data['abstract'] = get_abstract(article_el)

    data['author_list'] = get_authors(article_el)
    data['keyword_list'] = get_keywords(article_el)
    data['ref_list'] = get_references(article_el)

    return data


def parse_pubmed_xml(xml_name, out_json_name):

    print('reading ' + xml_name, flush=True)
    if xml_name.endswith('.gz'):
        data = gzip.open(xml_name).read().decode('utf-8')
    else:
        data = open(xml_name).read()

    print('parsing XML', flush=True)
    root = etree.XML(data.encode('utf-8'))
    article_list = root.getchildren()

    print('writing ' + out_json_name, flush=True)
    fout = open(out_json_name, 'w')
    narticle = 0
    for article in article_list:
        data = get_article_info(article)
        print(data['pmid'], json.dumps(data, sort_keys=True, ensure_ascii=False),
              sep='\t', file=fout, flush=True)
        narticle += 1
        if narticle % 1000 == 0:
            print('processed {} articles'.format(narticle), flush=True)
        # if narticle == 10:
        #     break

    fout.close()


if __name__ == '__main__':

    xml_name, out_json_name = sys.argv[1:3]
    parse_pubmed_xml(xml_name, out_json_name)
