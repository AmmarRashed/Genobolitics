try:
    from functools import lru_cache
except ImportError:
    from repoze.lru import lru_cache

import pyhgnc
import GEOparse
import numpy as np
from joblib import Parallel, delayed


@lru_cache(maxsize=None)
def lookup_gene(gene_symbol):
    try:
        query = lookup_gene.query
    except Exception as e:
        lookup_gene.query = query = pyhgnc.query()

    return query.hgnc(symbol=gene_symbol)


counter = 0

def get_genes_fold_changes(sample, label=None, sample_name=None, total=None):
    global counter
    genes_fold_changes = dict()
    for gene in np.unique(sample.index):
        hgnc = lookup_gene(gene)

        if hgnc:
            genes_fold_changes["HGNC:{}".format(hgnc[0].identifier)] = np.median(sample.loc[gene])
    counter += 1
    if counter % 25 == 0 or counter == total:
        print("{}:\t{} \t out of {}".format(counter, sample_name, total))

    if label is None:
        return genes_fold_changes
    return genes_fold_changes, label


def parse_database(geo_database, labels, geo=True, index_column="IDENTIFIER", prallel=True):
    global counter
    counter = 0
    if geo:
        df = GEOparse.get_GEO(geo=geo_database).table
    else:
        df = geo_database

    df = df.dropna().set_index(index_column)
    n = len(labels)
    X, y = list(range(n)), list(range(n))
    if prallel:
        data = Parallel(n_jobs=-1)(delayed(get_genes_fold_changes)(df[[sample]], label, sample, n) for sample, label in
                        labels.items())
        for i in range(n):
            X[i] = data[i][0]
            y[i] = data[i][1]
    else:
        i = 0
        for sample, label in labels.items():
            gene_fold_change, _ = get_genes_fold_changes(df[[sample]], label, sample, n)
            X[i] = gene_fold_change
            y[i] = label
            i += 1
    return X, y

