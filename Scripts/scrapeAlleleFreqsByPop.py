import urllib
from bs4 import BeautifulSoup
import pandas as pd
import time
import pdb
from collections import defaultdict, Counter
"""
page wise download
Beware: BeatuifulSoup uses best installedHTML parser
lxml failed to parse entire table. using python's parser now
"""
def extractData(row):
    cells = row.findAll("td")
    if len(cells) > 5:
        c1 = cells[1]
        alleleID = c1.find("a")['href'].split("=")[-1]
        try:
            alleleFreq = float(cells[5].text)
        except ValueError:
            print "Warning: allele freq not a float (assuming 0)", cells[3].text
            print row
            return alleleId, 0
        return alleleID, alleleFreq
    else:
        print "Warning, row to short", len(cells)
        print row
        return None, None

#alleleFreqDicts = []
with open("popListAFnet.pcl") as f: popList = cPickle.load(f)
"""http://www.allelefrequencies.net/hla6006a.asp?hla_locus_type=Classical&hla_population=2730"""
for pop in popList:
    print pop
    page = 1
    afsaccu = defaultdict(float)
    while True:

        url = 'http://www.allelefrequencies.net/hla6006a.asp?page=%s&hla_population=%s' % (page, pop)
        html = urllib.urlopen(url).read()
        soup = BeautifulSoup(html, "html.parser")

        tables = soup.find_all("table")
        entries = tables[4].findAll("tr")[1:]
        for row in entries:
            allele, af = extractData(row)
            if af>0:
                afsaccu[tuple(allele.split(":"))] = af
            ## accumulating allele freqs by HLA id-first 2 digits
            #if not allele is None and allele[0] in "ABC": 
            #    afsaccu[allele.split(":")[0].replace("*", "")] += af

        page += 1
        if not "hla6006a.asp?page=%s" % page in html:
            break
#    try:
#        afsaccu['PopSize'] = int(row.findAll("td")[7].text)
#    except ValueError, KeyError:
#        print 
#        afsaccu['PopSize'] = -1
    #afs = [(popId, af) for (popId, af) in afs if popId and not af is None]
        
    alleleFreqDicts.append(afsaccu)

df = pd.DataFrame(alleleFreqDicts, index=popList)

    
#df['Beduins'] = beduinData[0]

