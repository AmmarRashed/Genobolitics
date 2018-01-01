import pandas as pd, numpy as np
import GEOparse, pyhgnc, pickle, copy

def get_gene_fc_dict(gene_fc_dataframe_, sample_name, duplicate_handler=np.mean, hgnc=True, hgnc_path="", query=None):
    gene_fc_dataframe = copy.deepcopy(gene_fc_dataframe_)[[sample_name]]
    gene_fc_dict = dict()
    geo_hgnc_ = False
    try:
        geo_hgnc_dict = pickle.load(open("datasets/geo_to_hgnc_ids","rb"))
    except FileNotFoundError:
        geo_hgnc_ = True
        geo_hgnc_dict = dict()
    for row in gene_fc_dataframe.iterrows():
        gene_symbol = row[0]
        if hgnc:
            try:
                if geo_hgnc_:
                    geo_hgnc_dict[gene_symbol] = "HGNC:{}".format(query.hgnc(symbol=gene_symbol)[0].identifier)
                gene_symbol = geo_hgnc_dict[gene_symbol]
            except (KeyError, IndexError) as e:
                continue
        gene_fc = row[1]
        gene_fc_dict.setdefault(gene_symbol, [])
        gene_fc_dict[gene_symbol].append(gene_fc)
    for gene in gene_fc_dict:
        gene_fc_dict[gene] = duplicate_handler(gene_fc_dict[gene])
    if geo_hgnc_ and hgnc_path:
        with open("datasets/geo_to_hgnc_ids","wb") as f:
            pickle.dump(geo_hgnc_dict, f)
    print ("Dict for {} addedd successfully with {} samples".format(sample_name, len(gene_fc_dict)))
    return gene_fc_dict

def parse_database(geo_database_name, index_column="IDENTIFIER", hgnc_path="", labels={}):
    geo_df = GEOparse.get_GEO(geo=geo_database_name)\
                                .table\
                                .set_index(index_column).dropna()
    geo_df = geo_df[[col for col in geo_df.columns if "GSM" in col]]
    query = pyhgnc.query()
    X = list()
    y = list()
    for sample in geo_df.columns:
        if sample in labels:
            X.append(get_gene_fc_dict(geo_df, sample, query=query, hgnc_path=hgnc_path))
            y.append(labels[sample])
    return X,y
    
