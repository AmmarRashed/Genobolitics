import pickle

from metabolitics.preprocessing import MetaboliticsPipeline
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest

from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score
from imblearn.over_sampling import RandomOverSampler

import numpy as np, pandas as pd
from collections import defaultdict
from itertools import chain
from sklearn.model_selection import KFold

class GenoClassifier:
    def __init__(self, results_path, labels_path, diff=True, scale=False, select_features=False):
        results = pickle.load(open(results_path, 'rb'))
        labels = pickle.load(open(labels_path, 'rb'))
        self.binarize = lambda ls: [1 if l == 'unhealthy' else 0 for l in ls]

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

        if select_features:
            self.feature_selection_pipe = Pipeline([('select_k_best', SelectKBest(k=100))])
        else:
            self.feature_selection_pipe = None

        if scale:
            std_scalar = StandardScaler().fit(dataset, dataset.index)
            self.X, self.y = std_scalar.transform(dataset), dataset.index
        else:
            self.X, self.y = dataset, dataset.index

        self.X, self.y = np.array(self.X), np.array(self.y)

        self.models = [
              (RandomForestClassifier, {
                'max_depth': range(5, 20),
                'n_estimators': range(1, 8)
              }),

              (LogisticRegression, {
                'C': np.geomspace(3e-6, 1e2, num = 10),
                'max_iter': range(100, 1000 + 1, 1000)
              }),

              (SVC, {
                'C': np.geomspace(3e-6, 1e2, num = 10),
                'degree': range(1, 5),

                'max_iter': range(50000, 100000 + 1, 10000),
                'kernel': ['rbf', 'linear']
              }),

              (SGDClassifier, {
                'penalty': ['l1', 'l2', 'elasticnet'],
                'alpha': np.geomspace(1e-4, 1.0, num = 10),
                'max_iter': range(200, 1000 + 1, 100)
              }),
              (MLPClassifier, {
                'max_iter':range(500, 1001, 100),
                'activation':['relu', 'logistic','tanh', 'identity']
              })
            ]

    def _pipeline(self, X_train, X_test, y_train, y_test):
        try:
            feature_selection = self.feature_selection_pipe.fit(X_train, y_train)
            X_train_f, y_train_f = RandomOverSampler(random_state=42). \
                fit_sample(feature_selection.transform(X_train), y_train)

            X_test_f = feature_selection.transform(X_test)
        except:
            X_train_f, y_train_f, X_test_f= X_train, y_train, X_test

        cv_estimators = []
        for model, params in self.models:
            cv_model = GridSearchCV(model(random_state=42), params, n_jobs=-1).fit(X_train_f, y_train_f)
            cv_estimators.append(cv_model)

        f1_scores = []

        for estimator in cv_estimators:
            score = f1_score(self.binarize(estimator.predict(X_test_f)), self.binarize(y_test))
            f1_scores.append(score)

        best = cv_estimators[np.argmax(f1_scores)]
        best_estimator = best.best_estimator_
        print(best_estimator, best.best_score_)

        metrics = {'recall': recall_score, 'precision': precision_score, 'f1': f1_score,
                   'accuracy': accuracy_score}
        res = dict()
        y_predicted = self.binarize(best_estimator.predict(X_test_f))
        y_test_b = self.binarize(y_test)
        for metric, f in metrics.items():
            res[metric] = f(y_predicted, y_test_b)

        return res

    def eval_model(self, model, X_train, y_train, X_test, y_test):
        model.fit(X_train, y_train)
        pred_y = self.binarize(model.predict(X_test))
        metrics = {'recall': recall_score, 'precision': precision_score, 'f1': f1_score,
                   'accuracy': accuracy_score}
        res = dict()
        for m, f in metrics.items():
            res[m] = f(pred_y, self.binarize(y_test))
        return res


    def classify(self, k=10, model=None):
        metrics = {'recall': 0, 'precision': 0, 'f1': 0,
                   'accuracy': 0}
        kf = KFold(n_splits=k, random_state=42)
        for train_index, test_index in kf.split(self.X, self.y):
            X_train, X_test = self.X[train_index], self.X[test_index]
            y_train, y_test = self.y[train_index], self.y[test_index]
            if model is None:
                pred = self._pipeline(X_train, X_test, y_train, y_test)
            else:
                pred = self.eval_model(model, X_train, y_train, X_test, y_test)
            for m in pred:
                metrics[m] += pred[m]
        for m in metrics:
            metrics[m] /= k
        return metrics

