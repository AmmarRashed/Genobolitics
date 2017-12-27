import pickle
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

import numpy as np

X = pickle.load(open("results2.pickl", "rb"))
y = pickle.load(open("datasets/gene_data_y","rb"))

models = [(LogisticRegression, {'C': np.geomspace(1e-4, 1e4, num=10), 'max_iter': range(200, 1000 + 1, 100)}),
          (SVC, {'C': np.geomspace(1e-4, 1e4, num=10),
                 'max_iter': range(200, 1000 + 1, 100), 'kernel': ['rbf', 'poly', 'linear', 'sigmoid']}),
          (RandomForestClassifier, {'max_depth': range(10, 20), 'n_estimator': range(8)})]

results = []
for model, params in models:
    print(params)
    result = GridSearchCV(model(random_state=42), params).fit(X, y)
    results.append(result.best_score_)

best_model = models[np.argmax(results)].best_model_