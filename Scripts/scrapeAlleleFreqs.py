"Using BeatifulSoup html parsing to scrape systematically all "

import urllib
from bs4 import BeautifulSoup
import pandas as pd
import time
import pdb
from collections import defaultdict
import re

popIdDict = {}
popSizeDict = defaultdict(list)

def extractData(row):
    cells = row.findAll("td")
    if len(cells) > 5:
        c1 = cells[1]
        popId = c1.find("a")['href'].split("=")[-1]
        popIdDict[popId] = c1.text
        try:
            alleleFreq = float(cells[3].text)
            popSizeDict[popId].append(int(cells[5].text))
        except ValueError:
            print "Warning: allele freq not a float (assuming 0)", cells[3].text
            print row
            return popId, 0
        return popId, alleleFreq
    else:
        print "Warning, row to short", len(cells)
        print row
        return None, None

alleleFreqDicts = []
locuslist = ['A*01', 'A*02', 'A*03', 'A*11', 'A*23', 'A*24', 'A*25', 'A*26', 'A*29', 'A*30', 'A*31', 'A*32', 'A*33', 'A*34', 'A*36', 'A*43', 'A*66', 'A*68', 'A*69', 'A*74', 'A*80', 'B*07', 'B*08', 'B*13', 'B*14', 'B*15', 'B*18', 'B*27', 'B*35', 'B*37', 'B*38', 'B*39', 'B*40', 'B*41', 'B*42', 'B*44', 'B*45', 'B*46', 'B*47', 'B*48', 'B*49', 'B*50', 'B*51', 'B*52', 'B*53', 'B*54', 'B*55', 'B*56', 'B*57', 'B*58', 'B*59', 'B*67', 'B*73', 'B*78', 'B*81', 'B*82', 'B*83', 'C*01', 'C*02', 'C*03', 'C*04', 'C*05', 'C*06', 'C*07', 'C*08', 'C*12', 'C*14', 'C*15', 'C*16', 'C*17', 'C*18']
#locuslist=['B*52']

for locus in locuslist: 
    url = 'http://www.allelefrequencies.net/hla6002a.asp?all_name=%s' % locus
    html = urllib.urlopen(url).read()
    soup = BeautifulSoup(html, 'html.parser') ## if lxml is used, this truncates long tables!! (silently!!)
    
    tables = soup.find_all("table")
    try:
        checksum = int(re.findall("Allele reported (\d+)", tables[2].findAll("tr")[-1].text)[0])
    except IndexError, ValueError:
        checksum = None

    afs = [extractData(row) for row in tables[3].findAll("tr")[1:]]
    afs = [(popId, af) for (popId, af) in afs if popId and not af is None]
    alleleFreqDicts.append(dict(afs))
    print locus, len(afs), checksum
    time.sleep(1)

df = pd.DataFrame(alleleFreqDicts, index=[l.replace("*", "") for l in locuslist])

#df.to_csv("completeAlleleFrequencyTable2.csv")
