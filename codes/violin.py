# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 9:02:21 2022

@author: PI
"""

import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns

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

# Violin plot for y sample data distribution
fz = 24
fl = "Arial"
i = 2

y_violin = pd.DataFrame({
    "Data distribution of Cutoff": [""],
    "lambda cutoff(nm)": list(y)
})
ax = sns.violinplot(x="lambda cutoff(nm)", y="Data distribution of Cutoff", data=y_violin)

# Customize plot appearance
ax = plt.gca()
ax.spines['left'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['top'].set_linewidth(3)
ax.figure.set_size_inches(15, 10)
plt.grid(linewidth=0.3)
plt.yticks(fontproperties=fl, size=fz-2)
plt.xticks(fontproperties=fl, size=fz-2)
plt.xlabel('λcutoff(nm)', family=fl, fontsize=fz)
plt.ylabel('Frequency', family=fl, fontsize=fz)

# Feature selection for Top15 important features, you can  do the same to Top50 (X_filtered_50)
Top15_Descriptor = ['qed', 'SMR_VSA4', 'RingCount', 'SlogP_VSA6', 'SlogP_VSA4', 'fr_Ndealkylation2', 'SlogP_VSA11',
                    'EState_VSA7', 'EState_VSA4', 'NumAromaticCarbocycles', 'fr_NH2', 'fr_ether', 'PEOE_VSA1',
                    'EState_VSA6', 'HeavyAtomMolWt']
important_feature = [Top15_Descriptor[5*i], Top15_Descriptor[5*i+1], Top15_Descriptor[5*i+2], Top15_Descriptor[5*i+3],
                     Top15_Descriptor[5*i+4]]
X_important_feature = X_normalized[important_feature]

# Combine values for plotting
combined_values = X_important_feature[important_feature[0]].append(X_important_feature[important_feature[1]],
                                                                   ignore_index=True)
combined_values = combined_values.append(X_important_feature[important_feature[2]], ignore_index=True)
combined_values = combined_values.append(X_important_feature[important_feature[3]], ignore_index=True)
combined_values = combined_values.append(X_important_feature[important_feature[4]], ignore_index=True)

# Create a new DataFrame with combined values
new_df = pd.DataFrame({'Normalized Values': combined_values})

# Add a 'feature' column to the new DataFrame
new_df['feature'] = pd.concat([pd.Series([important_feature[0]] * len(X_important_feature)),
                               pd.Series([important_feature[1]] * len(X_important_feature)),
                               pd.Series([important_feature[2]] * len(X_important_feature)),
                               pd.Series([important_feature[3]] * len(X_important_feature)),
                               pd.Series([important_feature[4]] * len(X_important_feature))], ignore_index=True)

# Violin plot for normalized feature data
ax = sns.violinplot(x="feature", y="Normalized Values", data=new_df, dpi=300)

# Customize plot appearance
ax = plt.gca()
ax.spines['left'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['top'].set_linewidth(3)
sns.set_context(rc={'patch.linewidth': 12.0})
ax.figure.set_size_inches(12, 6)
plt.yticks(fontproperties=fl, size=fz-8)
plt.xticks(fontproperties=fl, size=fz-10)
plt.xlabel('Feature', family=fl, fontsize=fz-6)
plt.ylabel('Normalized Values', family=fl, fontsize=fz-6)
plt.savefig(r'..\figures\violin\FigS2.png')
# Show the plot
plt.show()
