import pickle

from Genobolitics import *
from metabolitics.preprocessing import MetaboliticsTransformer

genobolitics_transformer = MetaboliticsTransformer()
genobolitics_transformer.analyzer = Genobolitics("recon2")

gene_data = pickle.load(open("datasets/gene_data","rb"))

results = genobolitics_transformer.transform(gene_data)

pickle.dump(results, open("results1.pickl", "wb"))
pickle.dump(results, open("results1_backup.pickl", "wb"))

import pdb
pdb.set_trace()