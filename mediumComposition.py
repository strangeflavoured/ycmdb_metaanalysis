import re, numpy as np

#molar masses in g/mol
m_BO3=58.81
m_Ca=40.08
m_Cl=35.45
m_Co=58.93
m_Cu=63.55
m_Fe=55.84
m_I=126.90
m_K=39.10
m_Mg=24.31
m_Mn=54.94
m_Mo=96
m_NH4=18.04
m_NO3=62.01
m_Na=22.99
m_Ni=58.69
m_PO4=94.97
m_SO4=96.07
m_Se=78.97
m_Zn=65.4

#replace molar masses of hydrates
def replace_hydrate(comp,value,m_original):
	m_H2O=18.015
	if re.search(r"hydrate", comp, re.IGNORECASE):		
		if re.search(r"\s*dihydrate\s*", comp, re.IGNORECASE):
			value-=2*m_H2O/m_original
			comp=re.sub(r"\s*dihydrate\s*","", comp, re.IGNORECASE)
		elif re.search(r"\s*trihydrate\s*", comp, re.IGNORECASE):
			value-=3*m_H2O/m_original
			comp=re.sub(r"\s*trihydrate\s*","", comp, re.IGNORECASE)
		elif re.search(r"\s*tetrahydrate\s*", comp, re.IGNORECASE):
			value-=4*m_H2O/m_original
			comp=re.sub(r"\s*tetrahydrate\s*","", comp, re.IGNORECASE)
		elif re.search(r"\s*pentahydrate\s*", comp, re.IGNORECASE):
			value-=5*m_H2O/m_original
			comp=re.sub(r"\s*pentahydrate\s*","", comp, re.IGNORECASE)
		elif re.search(r"\s*hexahydrate\s*", comp, re.IGNORECASE):
			value-=6*m_H2O/m_original
			comp=re.sub(r"\s*hexahydrate\s*","", comp, re.IGNORECASE)
		elif re.search(r"\s*heptahydrate\s*", comp, re.IGNORECASE):
			value-=7*m_H2O/m_original
			comp=re.sub(r"\s*heptahydrate\s*","", comp, re.IGNORECASE)
		elif re.search(r"\shydrate\s*|\s*monohydrate\s*", comp, re.IGNORECASE):
			value-=m_H2O/m_original
			comp=re.sub(r"\s*hydrate\s*|\s*monohydrate\s*","", comp, re.IGNORECASE)
		else:
			print("Hydrate could not be identified: ", comp)

	return comp, value

#extract stochiometric factor of an element from molecular formula
#returns float or np.nan
def getStochiometry(element, formula):
	try:
		#search for digits following the element in formula
		regex=re.compile(f"(?<={element})[0-9]*")
		n=re.search(regex,formula).group(0)
		if not n:
			n=1.0
		else:
			n=float(n)
	#occurs if search has no result i.e. element is absent
	except AttributeError:
		n=0
	return n

#adds ion amount to total
#removes ion name from component name
#removes ion amount from orig. component
#returns total, adjusted component name and value
def getOneIon(total,total_name,comp,orig_comp,symbol,name,formula,value,mask,m_original,mass):
	if re.search(re.compile(name,re.IGNORECASE), comp):
		n=getStochiometry(symbol, formula)
		if n: 
			total[total_name][mask]+=value*n*mass/m_original
			value-=value*n*mass/m_original
			comp=re.sub(re.compile(f"\s*[A-z]*{name}[A-z]*\(*i*\)*", re.IGNORECASE),"", comp).strip()
	return total, comp, value

#caution: NH4, NO3 might give wrong results if organic compounds contain N
#also SO4; might also give wrong results with Se (but unprobable)
#returns updated total columns for input component
def getIons(m_original, medium, row, total):	
	comp=row["Component"]
	orig_comp=row["Component"]
	formula=row["Formula"]
	value=row["Value_g_l"]
	total_names=total.columns[1:]
	symbols=['B', 'Ca', 'Cl', 'Co', 'Cu', 'Fe', 'I', 'K', 'Mg', 'Mn', 'Mo',
		'N', 'N', 'Na', 'Ni', 'P', 'S', 'Se', 'Zn']
	names=["bor","calci","chlor","coba","cupr|copper","ferr|iron","iodi",
		"potass", "magne","manga","molyb","ammoni","nitra","sodi","nick",
		"phosphate", "sulfate","sele","zinc"]
	expr="|".join(names)
	#molar masses in g/mol
	masses=[m_BO3, m_Ca, m_Cl, m_Co, m_Cu, m_Fe, m_I, m_K, m_Mg, m_Mn, m_Mo,
		m_NH4, m_NO3, m_Na, m_Ni, m_PO4, m_SO4, m_Se, m_Zn]
	#check for expression in medium component
	#if re.search(re.compile(expr, re.IGNORECASE), comp):

	mask=total["Medium_uniqueID"]==medium
	#extract concentration for ions
	for i in range(len(masses)):
		total, comp, value=getOneIon(total,total_names[i],comp,orig_comp,symbols[i],
			names[i],formula,value,mask,m_original,masses[i])
		if value in [np.nan]:
			print(orig_comp,comp)

	orig=comp
	comp,value=replace_hydrate(comp,value,m_original)
	comp=cleanup(comp)
	
	if comp:
		row["Component"]=comp
		row["Value_g_l"]=value
	else:
		row=np.repeat(np.nan,len(row))	
	
	return row, total

def cleanup(comp):
	if comp.strip() in ["acid","ic", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa"]:
		comp=""
	else:
		comp=re.sub(re.compile(f"\s*([a-z]*hydro[A-z]*|ous)\s*", re.IGNORECASE),"", comp).strip()
	return comp