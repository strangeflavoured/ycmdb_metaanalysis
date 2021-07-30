#!/usr/bin/env python
#analysis routines for YCMDB database
#implemented: glucose uptake data

import config
from analysis import analysis, feature_matrix, pairplot

#sql connection
path=f"mysql+pymysql://{config.USER}:{config.PASSWORD}@{config.HOST}/{config.DB_NAME}"

#run analysis
data, strain_matrix, media_matrix=analysis(path,clean_unit=True,full=True)

features=feature_matrix(data, media_matrix, strain_matrix)

pairplot(features)


#if required: save data as csv
#with open("analysis.csv", "w+") as file:
#	file.write(data.to_csv())