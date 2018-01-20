from collections import defaultdict, OrderedDict
from itertools import chain, product
from functools import reduce

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelBinarizer
from sklearn.model_selection import GridSearchCV, cross_validate, StratifiedKFold
from sklearn.preprocessing import StandardScaler

import numpy as np, pandas as pd

from genobolitics import *

from metabolitics.preprocessing import MetaboliticsPipeline


def build_pipeline(p):
    pipeline, pipeline_params = [], OrderedDict()

    for model, model_params in p:
        name = model.__name__

        pipeline.append((name, model()))
        pipeline_params.update({'{}__{}'.format(name, param_name): values
                                for param_name, values in model_params.items()})

    return Pipeline(pipeline), pipeline_params


def build_pipelines(*args):
    return list(map(build_pipeline, product(*args)))


def preprocess_results(results, labels, use_diff_score, use_pathways, scale=True, use_one_hot=True):
    metabolitics_pipeline = []
    stage_1, scaled, binarized = None, None, None

    if use_diff_score:
        metabolitics_pipeline.append('reaction-diff')

    if use_pathways:
        metabolitics_pipeline.append('pathway_transformer')

    if metabolitics_pipeline:
        stage_1 = MetaboliticsPipeline(metabolitics_pipeline).fit_transform(results, labels)

    dataframe = get_dataframe(stage_1 or results, labels)

    if scale:
        scaled = StandardScaler().fit_transform(dataframe)

    if use_one_hot:
        binarized = LabelBinarizer().fit_transform(labels).ravel()

    return reduce(lambda x, y: y if y is not None else x, [dataframe, scaled]), \
           reduce(lambda x, y: y if y is not None else x, [labels, binarized])


def get_dataframe(X_dict, labels):
    samples = defaultdict(lambda: [])
    [
        samples[key].append(value) for key, value in
        chain(*map(lambda sample: sample.items(), X_dict))
    ]

    return pd.DataFrame(samples, index=labels)


def nested_cross_validation(X, y, pipelines, num_trials=10):
    metrics = ['f1', 'recall', 'precision', 'accuracy']
    trials = []

    for i in range(num_trials):
        cv_pipelines = []
        inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=i)
        outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=i)

        for pipeline, params in pipelines:
            cv_pipeline = GridSearchCV(pipeline, params, cv=inner_cv, n_jobs=-1, verbose=1).fit(X, y)
            cv_pipelines.append(cv_pipeline)

        best_pipeline = cv_pipelines[np.argmax([cv.best_score_ for cv in cv_pipelines])]
        cv = cross_validate(best_pipeline.best_estimator_,
                            X=X, y=y, cv=outer_cv, scoring=metrics,
                            return_train_score=False)

        trials.append((best_pipeline, cv))
        print("{} trial done\n{}".format(i + 1, '-' * 10))

    trials_scores = [scores for model, scores in trials]
    trials_means = map(lambda trial_scores: {key: value.mean()
                                             for key, value in trial_scores.items()}, trials_scores)

    return trials, pd.DataFrame(list(trials_means))
