import numpy as np
from functools import reduce, lru_cache

from metabolitics.analysis import MetaboliticsAnalysis
import GEOparse, pyhgnc


class FoldChange:
    def __init__(self, fold_change):
        self.fold_change = fold_change

    # OR IS MAX, MAX IS PLUS, THUS OR IS PLUS
    def __add__(self, other):
        return FoldChange(max(self.fold_change, other.fold_change))

    # AND IS MIN, MIN IS MINUS, THUS AND IS MINUS
    def __sub__(self, other):
        return FoldChange(min(self.fold_change, other.fold_change))


class Genobolitics(MetaboliticsAnalysis):
    def __init__(self, *args, **kwargs):
        super(Genobolitics, self).__init__(*args, **kwargs)
        self.model.solver = 'cplex'
        self.model.solver.configuration.timeout = 10 * 60

    def set_objective(self, measured_genes):
        self.clean_objective()
        for r in self.model.reactions:
            fold_change = 0
            
            try:
                fold_change = self.get_reaction_fold_change(r, measured_genes)
            except Exception as e:
                print('could not calculate fold changes for {}'.format(r))
            
            r.objective_coefficient = fold_change

    def get_reaction_fold_change(self, reaction, measured_genes):
        op = [('or', '+'), ('and', '-')]
        genes = [(gene, 'FoldChange({})'.format(self.get_gene_fold_change(gene, measured_genes))) for gene in
                 self.get_reaction_genes(reaction)]

        expr = reduce(lambda x, y: x.replace(*y), op + genes, reaction.gene_reaction_rule)

        return eval(expr).fold_change if expr != '' else 0

    def get_gene_fold_change(self, gene, measured_genes):
        return measured_genes[gene]

    def get_reaction_genes(self, reaction):
        return [g.id for g in reaction.genes]


@lru_cache(maxsize=None)
def lookup_gene(gene_symbol):
    try:
        query = lookup_gene.query
    except Exception as e:
        lookup_gene.query = query = pyhgnc.query()

    return query.hgnc(symbol=gene_symbol)


def get_genes_fold_changes(sample):
    genes_fold_changes = dict()
    for gene in np.unique(sample.index):
        hgnc = lookup_gene(gene)

        if hgnc:
            genes_fold_changes["HGNC:{}".format(hgnc[0].identifier)] = np.median(sample.loc[gene])

    return genes_fold_changes


def parse_database(geo_database_name, labels, index_column="IDENTIFIER"):
    geo_df = GEOparse.get_GEO(geo=geo_database_name).table.dropna().set_index(index_column)
    samples_columns = [column for column in geo_df.columns if column.startswith('GSM')]

    X, y = [], []
    for sample in labels.keys():
        gene_fold_change = get_genes_fold_changes(geo_df[[sample]])
        X.append(gene_fold_change)
        y.append(labels[sample])
        print("{} added with length {}".format(sample, len(gene_fold_change)))

    return X, y
