import pickle

from metabolitics.preprocessing import MetaboliticsPipeline

from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import  SelectKBest

from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score
from imblearn.over_sampling import RandomOverSampler

import numpy as np, pandas as pd
from collections import defaultdict
from itertools import chain

class GenoClassifier:
    def __init__(self, results_path, labels_path, diff=True):
        results = pickle.load(open(results_path, 'rb'))
        labels = pickle.load(open(labels_path, 'rb'))

        if diff:
            pipe = MetaboliticsPipeline(['reaction-diff',
                                         'pathway_transformer'])

            pre_processed_results = pipe.fit_transform(results, labels)
        else:
            pre_processed_results = results

        samples = defaultdict(lambda: [])
        [
            samples[key].append(value) for key, value in
            chain(*map(lambda sample: sample.items(), pre_processed_results))
        ]

        dataset = pd.DataFrame(samples, index=labels)

        std_scalar = StandardScaler().fit(dataset, dataset.index)
        self.feature_selection_pipe = Pipeline([('select_k_best', SelectKBest(k=100))])
        self.X, self.y = std_scalar.transform(dataset), dataset.index

        self.models = [
              (RandomForestClassifier, {
                'max_depth': range(5, 20),
                'n_estimators': range(1, 8)
              }),

              (LogisticRegression, {
                'C': np.geomspace(1e-4, 1e3, num = 10),
                'max_iter': range(100, 1000 + 1, 1000)
              }),

              (SVC, {
                'C': np.geomspace(1e-4, 1e2, num = 10),
                'degree': range(1, 5),

                'max_iter': range(50000, 100000 + 1, 10000),
                'kernel': ['rbf', 'linear']
              }),

              (SGDClassifier, {
                'penalty': ['l1', 'l2', 'elasticnet'],
                'alpha': np.geomspace(1e-4, 1.0, num = 10),
                'max_iter': range(200, 1000 + 1, 100)
              })
            ]


    def _pipeline(self, X_train, X_test, y_train, y_test):
        feature_selection = self.feature_selection_pipe.fit(X_train, y_train)
        X_train_f, y_train_f = RandomOverSampler(random_state=42). \
            fit_sample(feature_selection.transform(X_train), y_train)

        X_test_f = feature_selection.transform(X_test)
        cv_estimators = []
        for model, params in self.models:
            cv_model = GridSearchCV(model(random_state=42), params, n_jobs=-1).fit(X_train_f, y_train_f)
            cv_estimators.append(cv_model)

        f1_scores = []
        binarize = lambda ls: [1 if l == 'unhealthy' else 0 for l in ls]

        for estimator in cv_estimators:
            score = f1_score(binarize(estimator.predict(X_test_f)), binarize(y_test))
            f1_scores.append(score)

        best_estimator = cv_estimators[np.argmax(f1_scores)].best_estimator_

        metrics = {'recall': recall_score, 'precision': precision_score, 'f1': f1_score,
                   'accuracy': accuracy_score}
        res = dict()
        x_predicted = binarize(best_estimator.predict(X_test_f))
        y_test_b = binarize(y_test)
        for metric, f in metrics.items():
            res[metric] = f(x_predicted, y_test_b)

        return res

    def classify(self, k=10):
        metrics = {'recall': 0, 'precision': 0, 'f1': 0,
                   'accuracy': 0}
        for i in range(k):
            X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=0.10, stratify=self.y)
            pred = self._pipeline(X_train, X_test, y_train, y_test)
            for m in pred:
                metrics[m] += pred[m]
        for m in metrics:
            metrics[m] /= k
        return metrics
