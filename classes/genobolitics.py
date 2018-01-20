import math
from functools import reduce
from metabolitics.analysis import MetaboliticsAnalysis
import warnings


class FoldChange:
    def __init__(self, fold_change):
        self.fold_change = fold_change

    # OR IS MAX, MAX IS PLUS, THUS OR IS PLUS
    def __add__(self, other):
        return FoldChange(FoldChange.max_with_missing_values(self.fold_change, other.fold_change))

    # AND IS MIN, MIN IS MINUS, THUS AND IS MINUS
    def __sub__(self, other):
        return FoldChange(FoldChange.min_with_missing_values(self.fold_change, other.fold_change))

    @staticmethod
    def max_with_missing_values(*args):
        return max(FoldChange.replace_missing(*args, replace_with=-math.inf))

    @staticmethod
    def min_with_missing_values(*args):
        return min(FoldChange.replace_missing(*args, replace_with=math.inf))

    @staticmethod
    def replace_missing(*args, replace_with):
        replaced = [replace_with if o is None else o for o in args]
        if not all(args):
            warnings.warn('some operands are missing from logical expression!')
        return replaced


class Genobolitics(MetaboliticsAnalysis):
    def __init__(self, *args, **kwargs):
        super(Genobolitics, self).__init__(*args, **kwargs)
        self.model.solver = kwargs.get('solver', 'cplex')
        self.model.solver.configuration.timeout = kwargs.get('timeout', 10 * 60)

    def set_objective(self, measured_genes):
        warnings.filterwarnings("default")

        self.clean_objective()
        for r in self.model.reactions:
            fold_change = self.get_reaction_fold_change(r, measured_genes)

            if fold_change in {math.inf, -math.inf, None}:
                # missing operands (math.inf, -math.inf) or missing reaction rules (None)!
                r.objective_coefficient = 0.0
                warnings.warn('could not evaluate boolean expression, objective-coeff is set to ZERO!')
            else:
                r.objective_coefficient = fold_change

    def get_reaction_fold_change(self, reaction, measured_genes):
        op = [('or', '+'), ('and', '-')]
        genes_fold_changes = [(gene, 'FoldChange({})'.format(self.get_gene_fold_change(gene, measured_genes)))
                              for gene in self.get_reaction_genes(reaction)]

        expr = reduce(lambda x, y: x.replace(*y), op + genes_fold_changes, reaction.gene_reaction_rule)

        fold_change = None
        if expr:
            fold_change = eval(expr).fold_change
        else:
            # model returned the empty string as reaction rule
            warnings.warn('model returned the empty string as reaction rule!')

        return fold_change

    def get_gene_fold_change(self, gene, measured_genes):
        return measured_genes.get(gene, 'None')

    def get_reaction_genes(self, reaction):
        return [g.id for g in reaction.genes]
