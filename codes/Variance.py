# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:04:23 2022

@author: PI
"""

from sklearn.preprocessing import  MinMaxScaler
import pandas as pd


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

# Save variances to a CSV file
variances.to_csv(r'..\data\series_data.csv', index=True)
