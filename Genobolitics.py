import numpy as np

from metabolitics.analysis import MetaboliticsAnalysis

from functools import reduce


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
    def set_objective(self, measured_genes):
        self.clean_objective()
        for r in self.model.reactions:
            r.objective_coefficient = self.get_reaction_fold_change(r, measured_genes, duplicate_strategy=np.mean)

    def get_reaction_fold_change(self, reaction, measured_genes, duplicate_strategy):
        op = [('or', '+'), ('and', '-')]
        genes = [
            (gene, 'FoldChange({})'.format(self.get_gene_fold_change(gene, measured_genes, duplicate_strategy, nan=-1)))
            for gene in self.get_reaction_genes(reaction)]
        expr = reduce(lambda x, y: x.replace(*y), op + genes, reaction.gene_reaction_rule)
        if expr == "":
            return 0.0
        return eval(expr).fold_change

    def get_gene_fold_change(self, gene, measured_genes, duplicate_strategy=np.mean, nan=-1):
        try:
            return duplicate_strategy(measured_genes[gene])
        except KeyError:
            return nan

    def get_reaction_genes(self, reaction):
        return [g.id for g in list(reaction.genes)]



# measured_genes = pickle.load(open("datasets/gene_data","rb"))
# measured_genes_sample_1 = measured_genes[0]
# gb = Genobolitics("recon2")
# gb.set_objective(measured_genes_sample_1)