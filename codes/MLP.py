# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 13:12:14 2022

@author: PI
"""


import shap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import learning_curve

# Read X dataset from the Excel file
PI_cutoff = pd.read_excel(r"..\data\Descriptors.xlsx", sheet_name="descriptors_data")

# Extract features from the columns of the DataFrame
features = PI_cutoff.columns[1:]

# Extract X data
X = PI_cutoff.iloc[:, 1:]

# Check for missing values in the DataFrame
null_columns = X.columns[X.isnull().any()]
print("Columns with missing values:", null_columns)

# Convert X to float type
X = X.astype(float)

# Fill missing values with mean
for column in null_columns:
    mean_value = X[column].mean()
    X[column].fillna(mean_value, inplace=True)

# Recheck for missing values after filling
null_columns = X.columns[X.isnull().any()]
print("Columns with missing values after filling:", null_columns)

# Feature selection: Variance filtering
# Normalize data
scaler = MinMaxScaler()
X_normalized = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)

# Calculate variance for each column after normalization
variances = X_normalized.var()

# Select columns with variance greater than or equal to 0.05
selected_columns = variances[variances >= 0.05].index

# Keep the selected columns
X_filtered_50 = X[selected_columns]

# Assign X to the filtered dataset
X = X_filtered_50.values

# Read y dataset from the Excel file
PI_cutoffdata = pd.read_excel(r"..\data\CutoffData.xlsx", sheet_name="Sheet1")
y = PI_cutoffdata.iloc[:, 2]

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create a pipeline
pipeline = make_pipeline(StandardScaler(), MLPRegressor(random_state=42))

# Define hyperparameter candidate values
param_grid = {
    'mlpregressor__hidden_layer_sizes': [(10,), (20,), (30,), (40,), (10, 10),(10, 20),(10, 30),(10, 40), (10, 45),(20, 30), (30, 40), (40, 50)],
    'mlpregressor__activation': ['relu', 'tanh'],
    'mlpregressor__solver': ['sgd', 'adam'],
    'mlpregressor__alpha': [0.0001, 0.001, 0.01],
    'mlpregressor__learning_rate': ['constant', 'invscaling', 'adaptive'],
    'mlpregressor__learning_rate_init': [0.0001, 0.001, 0.01],
    'mlpregressor__max_iter': [1000, 2000, 5000, 10000],
    'mlpregressor__beta_1': [0.9, 0.95, 0.99],
    'mlpregressor__beta_2': [0.999, 0.9999]
}

# Use GridSearchCV for grid search
grid_search = GridSearchCV(pipeline, param_grid, cv=10, scoring='neg_mean_squared_error')
grid_search.fit(X, y)

best_model = grid_search.best_estimator_

# SHAP summary plot for feature importance
explainer = shap.KernelExplainer(best_model.predict, X_train)
shap_values = explainer.shap_values(X, nsamples=10000)

fig, ax = plt.subplots(dpi=300)
shap.summary_plot(shap_values, X, plot_type="bar", feature_names=selected_columns, color="c", max_display=15,
                  plot_size=(8, 8))
ax.spines['right'].set_visible(True)
ax.spines['top'].set_visible(True)
ax.spines['left'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['top'].set_linewidth(3)
plt.savefig(r'..\figures\Feature_importance.png', dpi=300)

# SHAP summary plot
fig, ax = plt.subplots(dpi=300)
shap.summary_plot(shap_values, X, feature_names=selected_columns, max_display=15, plot_size=(8, 5))
ax.spines['right'].set_visible(True)
ax.spines['top'].set_visible(True)
ax.spines['left'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['top'].set_linewidth(3)

ax.tick_params(axis='x', which='both', width=5)
ax.tick_params(axis='y', which='both', width=2)

x_ticklabels = ax.get_xticklabels()
for label in x_ticklabels:
    label.set_fontproperties('Arial')
    label.set_fontsize(14)

y_ticklabels = ax.get_yticklabels()
for label in y_ticklabels:
    label.set_fontproperties('Arial')
    label.set_fontsize(100)

plt.savefig(r'..\figures\shap_value.png', dpi=300)
plt.show()

# Plot learning curve
train_sizes, train_scores, val_scores = learning_curve(best_model, X_train, y_train, cv=5)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
val_scores_mean = np.mean(val_scores, axis=1)
val_scores_std = np.std