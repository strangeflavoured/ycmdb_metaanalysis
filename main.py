#!/usr/bin/env python
#analysis routines for YCMDB database
#implemented: glucose uptake data

from analysis import analysis, feature_matrix, pairplot, feature_pca

#run analysis
data, strain_matrix, media_matrix=analysis(clean_unit=True,full=False)

features=feature_matrix(data, media_matrix, strain_matrix)

#pairplot(features)

feature_pca(features)

#if required: save data as csv
#with open("analysis.csv", "w+") as file:
#	file.write(data.to_csv())