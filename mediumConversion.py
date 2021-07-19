#converts media information to a matrix, one column for each feature

import numpy as np
import pandas as pd
import re
from mediumComposition import getIons, getStochiometry

#to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None

#convert strings to integers
#undefined specifies how to deal with ValueErrors
#used for pubchem; molecular formula and unit after regex operation
#returns int, nan or none
def StringToInt(string, undefined="nan"):
	try:
		val=int(string)
	except ValueError:
		if undefined=="null":
			val=0
		elif undefined.lower() in ["nan", "na"]:
			val=float('NaN')
		else:
			val=None
	return val

#evaluates unit conversion formula
#undefined specifies how to treat NameErrors
#returns float, nan or none
def Eval(expression, undefined="NaN"):
	try:
		val=eval(expression)
	except NameError:
		if undefined in ["NaN", "nan", "NA", "na"]:
			val=float('NaN')
		elif undefined=="keep":
			val=expression
		elif undefined=="null":
			val=0
		else:
			val=None
	return val

#main function
#pulls media composition from database
def mediumConversion(con):
	# get list of medium IDs
	MediumIDs=pd.read_sql_query("select Medium_ID from RELATION_Medium", con)
	MediumIDs=MediumIDs["Medium_ID"].unique()

	# get medium components
	MediumTable=pd.read_sql_query("select * from META_MediumComposition", con)
	#replace na with np.nan to enable conversion to int
	MediumTable=MediumTable.fillna(np.nan)

	# replace ambiguities in unit characters
	MediumTable["Unit"]=MediumTable["Unit"].replace({"\xb5": "u"})	
	
	# convert pubchem to integer (or NA, this is fine, warning supressed)
	MediumTable["PubChem_orig"] = MediumTable["PubChem"] 
	pubchem=MediumTable["PubChem"].apply(lambda x: StringToInt(x))
	MediumTable["PubChem"]=pubchem

	# get lists of unique compounds, units and pubchems
	uniqueMediumCompounds = MediumTable["Component"].unique()
	uniqueUnits = MediumTable["Unit"].unique()
	uniquePubChemIDs = MediumTable["PubChem"].dropna().unique()
	uniquePubChemIDs =[str(int(x)) for x in uniquePubChemIDs]

	#pubchem Query molecular weight
	PubChemString=','.join(uniquePubChemIDs)
	pubChemUrl = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
		+PubChemString+
		"/property/MolecularWeight,MolecularFormula,IUPACName/CSV")
	uniquePubChemIDs = pd.read_csv(pubChemUrl)

	# extract number of carbons from formula	
	tmp = uniquePubChemIDs["MolecularFormula"]
	tmp = tmp.apply(lambda x: getStochiometry("C", x))
	uniquePubChemIDs["nC"]=tmp

	# original names in db of the unique compounds (also unify in MediumTable)
	uniquePubChemIDs["YCMNames"]=""
	for i, cid in enumerate(uniquePubChemIDs["CID"]):
		#mask where cid occurs in medium table
		idxCID = np.array([cid==pub for pub in MediumTable["PubChem"]])

		#get (first) name of components with cid
		component = MediumTable["Component"][idxCID].to_numpy()[0]

		#set names in medium table to selected name for cid
		MediumTable["Component"][MediumTable["PubChem"]==str(cid)] = component

		#set name for cid in uniquePubChemIDs table
		uniquePubChemIDs.at[i,"YCMNames"] = component

	# replace standard Media components	
	# MRL0000161  YNB
	# MRL0000162  YPD
	# MRL0000163  YPAD
	# MRL0000164  YPG
	# MRL0000165  SC medium
	# MRL0000166  YNB without ammonium sulfate
	# MRL0000167  YNB without ammonium sulfate and potassium
	# MRL0000171  SD medium
			  
	standardMedia =["MRL0000161", "MRL0000162", "MRL0000163", "MRL0000164",
		"MRL0000165", "MRL0000166", "MRL0000167", "MRL0000171"]

	#add standard medium composition for each medium which contains it
	#remove original rows which refer to standard medium
	for medium in standardMedia:
		# get the composition of the standard Medium		
		Composition = MediumTable[MediumTable["Medium_uniqueID"]==medium]

		# get all entries that contain this standard medium
		entries_ix = MediumTable.index[MediumTable["PubChem_orig"] == medium]
		for ix in entries_ix:
			RowsToAdd = Composition.copy()
			#associate the new rows with the current medium
			RowsToAdd["Medium_uniqueID"] = MediumTable["Medium_uniqueID"].iloc[ix]
			#add the new rows
			MediumTable = pd.concat([MediumTable, RowsToAdd])
			
		# delete initial rows
		MediumTable = MediumTable.drop(entries_ix)

	# unify names for peptone, yeast extract
	mask=MediumTable["PubChem_orig"]=="MRL0000170"
	MediumTable["Component"][mask] = "Peptone"
	peptone=MediumTable[mask]
	MediumTable=MediumTable[~mask]

	mask=MediumTable["PubChem_orig"]=="MRL0000169"
	MediumTable["Component"][mask] = "Yeast extract"
	yeast_extract=MediumTable[mask]
	MediumTable=MediumTable[~mask]

	# convert all units to g/l
	arr=np.repeat("",len(uniqueUnits))
	conversion = pd.DataFrame({"formula": arr}, index=uniqueUnits)

	conversion.loc["bool","formula"] =  ""
	conversion.loc["mM","formula"] =  "*1e-3 * molar_mass_g_mol"
	conversion.loc["uM","formula"] = "*1e-6 * molar_mass_g_mol"
	conversion.loc["\xb5M","formula"] = "*1e-6 * molar_mass_g_mol"
	conversion.loc["%","formula"] = "*10"   # assuming (w/v)
	conversion.loc["% (w/v)","formula"] = "*10"
	conversion.loc["g/l","formula"] = "*1"
	conversion.loc["mg/l","formula"] = "*1e-3"
	#conversion.loc["ml","formula"] = "*float('NaN')"
	conversion.loc["g/L","formula"] = "*1"
	conversion.loc["mg/L","formula"] = "*1e-3"
	conversion.loc["\xb5g/L","formula"] = "*1e-6"
	conversion.loc["ug/L","formula"] = "*1e-6"
	conversion.loc["M","formula"] = "* molar_mass_g_mol"
	conversion.loc["mg/ml","formula"] = "*1"
	conversion.loc["g/mol carbon","formula"] = "* mol_total_carbon_per_litre"
	conversion.loc["mg/mol carbon","formula"] = "*1e-3 * mol_total_carbon_per_litre"
	# density must de known, not available via pubchem easily 
	# *10e-3 * mol_total_carbon_per_litre"
	#conversion.loc["ml/mol carbon","formula"] = "*float('NaN')"
	#conversion.loc["% (v/v)","formula"] = "*float('NaN')" # mixed compound
	conversion.loc["ug/l","formula"] = "*1e-6"
	conversion.loc["g/1ml aqua dest.","formula"] = "*1e3"

	conversion.loc["NA","formula"] = "*float('NaN')"
	conversion.loc["","formula"] = "*float('NaN')"

	conversion.loc["umol/L","formula"] = "*1e-6 * molar_mass_g_mol"
	conversion[conversion.isna()] = "*NA"""

	# substitute molar_mass_g_mol
	# components w/o unit are dropped here
	MediumTable=MediumTable[~MediumTable["Unit"].isna()]

	#create and conversion formula for each component
	#formula is column of medium table
	conv=[]
	z=zip(MediumTable["Value"], conversion.loc[MediumTable["Unit"],"formula"])
	for i,j in z:
		conv.append(str(i)+str(j))
	MediumTable["conversion"] = conv
	#substitute molar masses in medium table conversion
	for indx, conv in enumerate(MediumTable["conversion"]):
		#find pubchem of row in pubchem table
		mask=uniquePubChemIDs["CID"]==MediumTable["PubChem"].iloc[indx]
		accessIndex=uniquePubChemIDs.where(mask).dropna().index
		if not accessIndex.empty:
			val=str(uniquePubChemIDs["MolecularWeight"][accessIndex].iloc[0])
		else:
			#component has no pubchem
			val="float('nan')"
		MediumTable["conversion"].iloc[indx]=conv.replace("molar_mass_g_mol", val)
	# deal with remaining mol_carbon units (number of carbons required)

	# find all entries with the unit 1/mol carbon, replace with "1"
	ixU = MediumTable["Unit"].apply(lambda x: re.sub(".*(/mol carbon).*","1", x))
	# find all media that contain compound with the unit
	ixM = MediumTable["Medium_uniqueID"][ixU=="1"].unique()
	# mask carbon sources: ethanol and glucose
	ixC = np.array(MediumTable["PubChem"].isin([702, 5793, 6255, 6036]))
  
	for MediumId in ixM:
		# find values for ethanol or glucose
		mask=np.array(MediumTable["Medium_uniqueID"]==MediumId) & ixC
		cSource = MediumTable[mask]

		#divide by molar weight, multiply wth number of carbons
		g_l = np.array(cSource["conversion"].apply(lambda x: Eval(x)))
		ixP = uniquePubChemIDs["CID"].isin(cSource["PubChem"])
		molWeight=np.array(uniquePubChemIDs["MolecularWeight"][ixP])
		nC=np.array(uniquePubChemIDs["nC"][ixP])
		molC = str(np.nansum(g_l/molWeight * nC))

		#substitute value for mol_total_carbon_per_litre
		mask=MediumTable["Medium_uniqueID"]==MediumId
		expr='mol_total_carbon_per_litre'
		conv=MediumTable["conversion"].apply(lambda x: x.replace(expr, molC))
		MediumTable["conversion"][mask] = conv
	
	# evaluate conversion formula
	conv=MediumTable["conversion"].apply(lambda x: Eval(x))
	MediumTable["Value_g_l"] =  conv
		
	# sort media into long list
	columns=["Medium_uniqueID","Value_g_l","Component","PubChem"]
	outmatrix = MediumTable[columns]
	# remove entries without ID
	outmatrix = outmatrix.dropna(subset=["PubChem"])
	pubChem=outmatrix[["Component","PubChem"]]
	molarMasses=uniquePubChemIDs[["CID","MolecularWeight"]]
	molarMasses=molarMasses.rename(columns={"CID": "PubChem"})
	molecularWeight=pd.merge(molarMasses, pubChem, how="left", on="PubChem")
	val=molecularWeight["MolecularWeight"].values
	molecularWeight=pd.Series(val, index=molecularWeight["Component"])
	molecularWeight=molecularWeight.to_dict()

	outmatrix["Formula"]=np.nan
	for i, row in outmatrix.iterrows():
		mask=uniquePubChemIDs["CID"]==row["PubChem"]
		formula=uniquePubChemIDs["MolecularFormula"][mask].to_numpy()[0]
		outmatrix.at[i,"Formula"]=formula

	#reindex
	outmatrix=outmatrix.reset_index()

	#sum ions	
	total_quantities=['Medium_uniqueID','total BO3', 'total Cl',
		'total Co', 'total Cu', 'total Fe', "total I",'total K', 'total Mg', 'total Mn',
		'total Mo', 'total NH4', 'total NO3', 'total Na','total Ca',
		'total Ni', 'total PO4', 'total SO4', 'total Se', 'total Zn']
	total=pd.DataFrame(columns=total_quantities)
	total["Medium_uniqueID"]=outmatrix["Medium_uniqueID"].unique()
	total=total.fillna(0)
	for medium in outmatrix["Medium_uniqueID"].unique():
		medium_mask=outmatrix["Medium_uniqueID"]==medium
		for i, row in outmatrix.iterrows():
			if medium_mask[i]==True:
				row, total =getIons(molecularWeight[row["Component"]], medium, row, total)
				outmatrix.iloc[i]=row
	outmatrix=outmatrix.dropna(subset=['Component'])	

	with open("outmatrix.csv", "w+") as file:
		file.write(outmatrix.to_csv())
	# remove Pubchem, Formula
	outmatrix = outmatrix[["Medium_uniqueID","Value_g_l","Component"]]
	for i, row in total.iterrows():
		for j in range(1,len(total.columns)):
			new_row={"Medium_uniqueID":row["Medium_uniqueID"], "Value_g_l":row[j], "Component":total.columns[j]}
			outmatrix=outmatrix.append(new_row, ignore_index=True)

	for ID in peptone["Medium_uniqueID"]:
		new_row={"Medium_uniqueID":ID, "Value_g_l":1, "Component":"Peptone"}
		outmatrix=outmatrix.append(new_row, ignore_index=True)
	for ID in yeast_extract["Medium_uniqueID"]:
		new_row={"Medium_uniqueID":ID, "Value_g_l":1, "Component":"Yeast Extract"}
		outmatrix=outmatrix.append(new_row, ignore_index=True)
	# spread stacked values to matrix
	#adds new column for each value of Component, unit is g/l
	outmatrix_wide = outmatrix.pivot_table(index="Medium_uniqueID", columns="Component", values="Value_g_l")
	outmatrix_wide=outmatrix_wide.fillna(0)

	return outmatrix_wide