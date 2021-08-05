import numpy as np
import pandas as pd

from mediumComposition import getIons
from mediumConversion import mediumConversion


def testIonCount():
	print("TEST ION COUNT")
	total_quantities=['Medium_uniqueID','total BO3', 'total Ca','total Cl',
		'total Co', 'total Cu', 'total Fe', "total I",'total K', 'total Mg', 'total Mn',
		'total Mo', 'total NH4', 'total NO3', 'total Na',
		'total Ni', 'total PO4', 'total SO4', 'total Se', 'total Zn']
	total=pd.DataFrame(columns=total_quantities)
	name="EXP0001234"
	total["Medium_uniqueID"]=[name]
	total=total.fillna(0)
	df=pd.DataFrame({"Component":["Ferrous ammonium sulfate hexahydrate"],"Formula":["FeH20N2O14S2"],"Value_g_l":[2]})
	row=df.iloc[0]
	molecular_weight=392.1
	row, total=getIons(molecular_weight, name, row, total)
	#expected in total: Fe 0.28483, NH4 0.18403, SO4 0.98006
	#row empty
	print(row)
	print("expected: empty")
	print(total.loc[:,(total != 0).any(axis=0)])
	print("expected: Fe 0.284825,  NH4 0.184035,  SO4 0.980056\n")

if __name__=="__main__":
	testIonCount()