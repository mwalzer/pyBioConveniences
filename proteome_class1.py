#!/usr/bin/env python
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
import numpy as np
import pandas as pd

import math
import multiprocessing
import pickle
import sys
import os
import time
import random
from collections import defaultdict

from Fred2.Core import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.Core import Allele

import gc
import argparse


__author__ = 'walzer'
VERSION = "0.0"


def random_peptide_sequence(length=9):
    return ''.join(random.choice("ACDEFGHIKLMNPQRSTVWY") for i in range(length))


def slice_mers(biopy_seqrec, windowsize=9, exclude_ambiguous=True):
    res = defaultdict(dict)
    agg = defaultdict(int)
    for item in biopy_seqrec:
        s = defaultdict(int)
        try:
            seq = str(item.seq)
            for i, v in enumerate(seq):
                if len(seq)-i < windowsize:
                    break
                tmp = seq[i:(i + windowsize)]
                if exclude_ambiguous:
                    if 'B' not in tmp and 'O' not in tmp and 'U' not in tmp and 'X' not in tmp and 'Z' not in tmp:
                        s[tmp] += 1
                        agg[tmp] += 1
                    else:
                        s[tmp] += 1
                        agg[tmp] += 1
        except:
            print 'error with item', item.name
        res[item.name]['len'] = len(item)
        for ntup, nval in enumerate(np.bincount(s.values())):
            res[item.name][str(ntup)+'tupel'] = nval
    for ntup, nval in enumerate(np.bincount(agg.values())):
        res["agglomerated"][str(ntup)+'tupel'] = nval
        gc.collect()
    return res, agg.keys()

    #k = slice_mers(records, 9, True, True)
    #sum(k["agglomerated"].values())


def toplevel_predictor(x):
    predictor = EpitopePredictorFactory("netMHC", version="3.4")
    return predictor.predict(x)


def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-t', '--threads', dest="threads", help='number of threads to use in parallel.', default=1, type=int)
    parser.add_argument('-l', '--length', dest="length", help='size of the ligands', default=9, type=int)
    parser.add_argument('-out', dest="out", help="<Required> full path to the output file", required=True, type=str)
    parser.add_argument('-f', '--fasta', dest="fasta", help='the fasta to be considered. Mutually exclusive to -w', type=str)
    parser.add_argument('-w', '--window', dest="window", help='the number of windows to be be considered (and created by random). Mutually exclusive to -f', type=int)
    parser.add_argument('-c', '--chunksize', dest="chunksize", help='max number of peptides subjected to prediction per thread at once', default=10000, type=int)

    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if options.fasta and options.window:
        parser.print_help()
        sys.exit(1)

    #http://stackoverflow.com/questions/5549190/is-shared-readonly-data-copied-to-different-processes-for-python-multiprocessing/5550156#5550156
    peps = list()
    #http://stackoverflow.com/questions/22487296/multiprocessing-in-python-sharing-large-object-e-g-pandas-dataframe-between

    #fastaname = "/share/usr/walzer/immuno-tools/dbs/swissprotHUMANwoi_130927.fasta"
    if options.fasta:
        with open(options.fasta, "rU") as handle:
            records = [record for record in SeqIO.parse(handle, "fasta")]
            #agl, peps = slice_mers(records, 9, True)
            peps = [Peptide(x) for x in slice_mers(records, 9, True)[1]]

    if options.window:
        for i in range(0, options.window):
            peps.append(Peptide(random_peptide_sequence(length=options.length)))

    gc.collect()

    import time
    start_time = time.time()

    pool = multiprocessing.Pool(processes=options.threads)
    chunks = [peps[x:x + options.chunksize] for x in xrange(0, len(peps), options.chunksize)]
    results = pool.map(toplevel_predictor, chunks)
    #results = [pool.apply(toplevel_predictor, args=(x,)) for x in chunks]

    #merge the dataframes does not work like supposed to! with same method DFs
    # result = results[0]
    # for x in results[1:]:
    #     result.merge_results(x)

    result = pd.concat(results)
    result.to_csv(options.out)

    tx = (time.time() - start_time)
    print("--- %s seconds ~ %s chunks---" % (tx, len(chunks)))


if __name__ == '__main__':
    __main__()