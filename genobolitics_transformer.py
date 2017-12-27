import pickle

from Genobolitics import *
from sklearn_utils.preprocessing import FoldChangeScaler
from sklearn.pipeline import Pipeline
from metabolitics.preprocessing import MetaboliticsTransformer, MetaboliticsPipeline

genobolitics_transformer = MetaboliticsTransformer()
genobolitics_transformer.analyzer = Genobolitics("recon2")

X = pickle.load(open("datasets/gene_data_X","rb"))
y = pickle.load(open("datasets/gene_data_y","rb"))
print(len(X),len(y))
pipe = Pipeline([
    ('scaling', FoldChangeScaler(reference_label='healthy')),
    ('LP_FVA', genobolitics_transformer),
    ('Metabolotics', MetaboliticsPipeline(
        ['reaction-diff',
         'feature-selection',
         'pathway_transformer']
    ))
])


results = pipe.fit_transform(X, y)
pickle.dump(results, open("results2.pickl", "wb"))
pickle.dump(results, open("results2_backup.pickl", "wb"))

import pdb
pdb.set_trace()
