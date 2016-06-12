#!/usr/bin/env python
__author__ = 'walzer'
import sys
import os
from datetime import datetime
import logging
import argparse
import pyopenms as oms
from os.path import join
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO import FileReader

VERSION = "0.3"

def categorize(score):
    if score >= 0.6384377847127609:
        return "strong"
    elif score >= 0.4256251898085073:
        return "weak"
    else:
        return None

def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-in',  dest="inf", help='<Required> full path to the input file', required=True)
    parser.add_argument('-out', dest="out", help="<Required> full path to the output file", required=True)
    parser.add_argument('-allele', dest="allele", help="<Required> full path to an allele file", required=True)

    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if not (options.inf or options.out):
        parser.print_help()
        sys.exit(1)

    target_alleles_set = set()
    if not (options.allele):
        #Fred2.FileReader.read_lines is broken
        #alleles = FileReader.read_lines(options.allele, type=Allele)
        with open(options.allele, 'r') as handle:
            for line in handle:
                target_alleles_set.add(Allele(line.strip().upper()))

    if not target_alleles_set:
        parser.print_help()
        sys.exit(1)
        
    ttn = EpitopePredictorFactory('netmhc')

    pros = list()
    peps = list()
    f = oms.IdXMLFile()
    f.load(options.inf, pros, peps)

    pepstr = set()
    for pep in peps:
        for h in pep.getHits():
            if "decoy" not in h.getMetaValue("target_decoy"):
                if 7 < len(h.getSequence().toUnmodifiedString()) < 12:
                    pepstr.add(h.getSequence().toUnmodifiedString())

    es = [Peptide(x) for x in pepstr]

    try:
        preds_n = ttn.predict(es, alleles=target_alleles_net_a)
    except Exception as e:
        print "something went wrong with the netMHC prediction", options.inf, "what:",  str(e)
        sys.exit(1)

    preds = dict()    
    for index, row in preds_n.iterrows():
        score = row.max()
        allele = str(row.idxmax())
        categ = categorize(score)
        seq = row.name[0].tostring()
        if categ:
            preds[seq] = (allele,categ,score)

    npeps = list()
    for pep in peps:
        hits = pep.getHits()
        nhits = list()
        for h in hits:
            if h.getSequence().toUnmodifiedString() in preds:
                x = preds[h.getSequence().toUnmodifiedString()]
                h.setMetaValue('binder',x[0])
                h.setMetaValue(str(x[1]),x[2])
                nhits.append(h)
            else:
                nhits.append(h)
        pep.setHits(nhits)                    

    f.store(options.out, pros, peps)


if __name__ == '__main__':
    __main__()