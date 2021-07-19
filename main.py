#!/usr/bin/env python
#analysis routines for YCMDB database
#implemented: glucose uptake data

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pymysql, csv, re
from sqlalchemy import create_engine
import pickle, math
import config

#imports from other files
from mediumConversion import mediumConversion

#pull and process data from database
#returns glucose uptake, medium and strain data
def analysis(clean_unit=False, full=False):

	#create mysql connection
	path=f"mysql+pymysql://{config.USER}:{config.PASSWORD}@{config.HOST}/{config.DB_NAME}"
	
	with create_engine(path).connect() as connection:
		#disable ONLY_FULL_GROUP_BY
		connection.execute("""set sql_mode=(select replace(@@sql_mode,
			'ONLY_FULL_GROUP_BY', ''))""")

		#get time-independent glucose uptake data
		query=("select Author_Year, Numerical_Value, Unit, Strain_ID, "
			"Medium_ID, Growth_Rate, Growth_Rate_Type, Growth_Phase, "
			"Synchronisation, Culture, Temperature, Aeration, Compartment, "
			"uniqueID from MET_Uptake where Compound='Glucose' and Time_Unit "
			"is null")
		data=pd.read_sql(query, connection)

		#get publication data from each publication in the data set
		conditions="' or Author_Year='".join(data["Author_Year"].unique())
		query=("select Author_Year, General_CellDensity, General_DryWeight "
			f"from RELATION_Publication where (Author_Year='{conditions}')")
		meta=pd.read_sql(query, connection)
		
		#get cell density for each publication
		authors=meta["Author_Year"][meta["General_CellDensity"]!=0].unique()
		conditions="' or Author_Year='".join(authors)
		query=("select Author_Year, Numerical_Value, Unit, Time, Time_Unit "
			f"from General_CellDensity where Author_Year='{conditions}'  "
			"and Time_Unit is null")
		meta_CD=pd.read_sql(query, connection)
		
		#get dry weight data for each publication
		authors=meta["Author_Year"][meta["General_DryWeight"]!=0].unique()
		conditions="' or Author_Year='".join(authors)
		query=("select Author_Year, Numerical_Value, Unit, Medium_ID, "
			"Growth_Rate from General_DryWeight where (Author_Year='"
			f"{conditions}') and Time_Unit is null")
		meta_DW=pd.read_sql(query, connection)
		
		#add dry weight data columns to data
		data["DryWeight"]=np.repeat(np.nan, data.shape[0])
		data["DW_Unit"]=np.repeat("", data.shape[0])


		#add DW data by mapping to growth rate and medium
		for i, row in data.iterrows():
			#map data
			GR_match = meta_DW["Growth_Rate"]==row["Growth_Rate"]
			Med_Match = meta_DW["Medium_ID"]==row["Medium_ID"]
			mask=np.all([GR_match, Med_Match], axis=0)
			DW=meta_DW[["Numerical_Value", "Unit"]][mask]
			try:
				#set values
				data.at[i,"DryWeight"]=DW["Numerical_Value"].to_numpy()[0]
				unit=DW["Unit"].to_numpy()[0]
				if unit in ["g/l", "mg/ml", "mg/mL"]:
					unit="g/L"
				data.at[i,"DW_Unit"]=unit
			except IndexError:
				continue

		#unit conversion for glucose uptake
		data=data.apply(lambda x: glucose_uptake_unit(x), axis="columns")

		#if clean unit mode, keep only data converted to g/gDWh
		if clean_unit:
			data=data[data["Unit"]=="g/(gDW h)"]

		#remove any points without numerical values
		data=data[data["Numerical_Value"].notnull()]	
		

		#slow, creates strain and media matrix
		if full:
			#create composition matrix for each medium contained in data
			#uses medium uniqueID for mapping

			#get medium relation for each medium
			media=data["Medium_ID"].unique()
			media="' or Medium_ID='".join(media)
			query=("select Medium_ID, uniqueID, pH from RELATION_Medium "
				f"where Medium_ID='{media}'")
			#use %% instead of % to escape % as formatter
			query=re.sub("%","%%", query)
			media=pd.read_sql(query,connection)
			
			#do some conversions on medium composition (c.f. mediumConversion.py)
			#return contains each medium's composition in single row
			media_matrix=mediumConversion(connection)
			print(media_matrix.columns)
			#restore uID from index to coĺumn
			media_matrix['uniqueID'] = media_matrix.index

			#merge medium relation and composition
			#only keeps media present in data
			#contains medID and uID, pH and composition
			media_matrix=pd.merge(media, media_matrix, how="left", on="uniqueID")

			#drop columns which contain no value
			media_matrix=media_matrix.dropna(axis="columns", how="all")

			#create strain info matrix for each strain in data
			#uses strain ID for mapping

			#get strain relation for strains in data
			strains=data["Strain_ID"].unique()
			strains="' or Strain_ID='".join(strains)
			query=f"select * from RELATION_Strain where Strain_ID='{strains}'"
			strains=pd.read_sql(query,connection)

			#only keep relevant columns (removes table counts)
			to_keep=['Strain_ID','Closest_Ancestor', 'Ploidity',
				'Mating_Type','Mutations', 'Mutations.1', 'Mutations.2', 
				'Mutations.3', 'Mutations.4', 'Mutations.5', 'Mutations.6',
				'Mutations.7', 'Mutations.8', 'Mutations.9', 'Mutations.10',
				'Mutations.11', 'Mutations.12']
			strains=strains[to_keep]

			#get a list of all occuring mutations
			mutations=['Mutations', 'Mutations.1', 'Mutations.2',
				'Mutations.3', 'Mutations.4', 'Mutations.5', 'Mutations.6',
				'Mutations.7', 'Mutations.8', 'Mutations.9', 'Mutations.10',
				'Mutations.11', 'Mutations.12']
			mutations=strains[mutations].values.flatten()
			mutations=list(np.unique(mutations[mutations != np.array(None)]))			

			#create strain matrix with all mutations
			columns=['Strain_ID','Closest_Ancestor','Ploidity','Mating_Type']
			strain_matrix=pd.DataFrame(columns=columns+mutations)
			strain_matrix[columns]=strains[columns]
			
			#check whether a mutation is present in a strain
			columns=strain_matrix.columns
			for i, row in strain_matrix.iterrows():
				#first mutation is in 4th column
				for j in range(4,len(row)):
					mask=strains["Strain_ID"]==row[0]
					if columns[j] in strains[mask].to_numpy():
						strain_matrix.iat[i,j]=True
					else:
						strain_matrix.iat[i,j]=False
				
			#drop columns which contain only the same value
			unique = strain_matrix.apply(pd.Series.nunique)
			where_same=unique[unique == 1].index
			strain_matrix=strain_matrix.drop(where_same, axis="columns")

			#pickle strain and media matrix as meta_matrix.pkl
			pickle.dump({"media_matrix": media_matrix, "strain_matrix":
				strain_matrix}, open("meta_matrix.pkl", "wb+"))

		#load pickled matrices if necessary
		else:
			try:
				meta=pickle.load(open("meta_matrix.pkl", "rb"))	
				strain_matrix=meta["strain_matrix"]
				media_matrix=meta["media_matrix"]
			except FileNotFoundError:
				print("No strain or media data found")
				strain_matrix=None
				media_matrix=None

	#return relevant data frames
	return data, strain_matrix, media_matrix

#converts glucose uptake units to g/gDWh if possible
#requires DryWeight to convert from /L to /gDW
#might be possible with cell density data, not implemented
#apply to each row individually, set returned obj as row
def glucose_uptake_unit(row, molar_mass=180.156):
	unit=row["Unit"]

	#mol->g
	if unit in ["mmol/(gDW h)","mmol/(g h)"]:
		row["Numerical_Value"]=row["Numerical_Value"]*1e-3*molar_mass
		row["Unit"]="g/(gDW h)"
	elif unit in ["umol/(gDW h)","umol/(g h)"]:
		row["Numerical_Value"]=row["Numerical_Value"]*1e-6*molar_mass
		row["Unit"]="g/(gDW h)"

	#decimal prefixes
	elif unit in ["mg/(mgDW h)", "g/(g h)"]:
		row["Unit"]="g/(gDW h)"

	#/L->/gDW
	elif unit in ["ug/ml","ug/mL"]:
		row["Numerical_Value"]=row["Numerical_Value"]*1e-3/row["DryWeight"]
		row["Unit"]="g/gDW"
	elif unit in ["mg/ml","mg/mL"]:
		row["Numerical_Value"]=row["Numerical_Value"]/row["DryWeight"]
		row["Unit"]="g/gDW"
	elif unit in ["umol/(ml s)","umol/(mL s)"]:
		row["Numerical_Value"]=row["Numerical_Value"]*1e-3*60**2/row["DryWeight"]
		row["Unit"]="g/(gDW h)"

	return row

if __name__ == '__main__':

	#run analysis
	data, strain_matrix, media_matrix=analysis(clean_unit=True,full=True)

	features=pd.merge(data, media_matrix, how="left", on="Medium_ID")
	features=pd.merge(features, strain_matrix, how="left", on="Strain_ID")
	features=features.drop(["Author_Year","Unit","Strain_ID","Medium_ID",
		"Growth_Rate","Growth_Rate_Type","Synchronisation","Compartment",
		"uniqueID_x","uniqueID_y","DryWeight","DW_Unit","Closest_Ancestor"],axis=1)
	for col in features.columns:
		if len(features[col].unique())<2:
			features=features.drop(col,axis=1)

	with open("features.csv", "w+") as file:
		file.write(features.to_csv())

	columns=list(features.columns)
	columns.remove("Numerical_Value")
	#columns.remove("Author_Year")
	sns.set_theme()
	for i in range(math.ceil(len(columns)/5)):
		cols=columns[5*i:5*(i+1)]
		sns.pairplot(features[cols+["Numerical_Value"]], y_vars=["Numerical_Value"],x_vars=cols)
		plt.savefig(f"images/pairplots_{i}.png")
		plt.close()

	#if required: save data as csv
	#with open("analysis.csv", "w+") as file:
	#	file.write(data.to_csv())
