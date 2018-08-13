import itertools
def aminotuple(repeat_val):
	x=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
	L=[]
	for p in itertools.product(x,repeat=repeat_val):
		L.append(p)

	return L

def aminodict():
	x={'ALA':['N','CA','C','O','CB'], 'VAL':['N','CA','C','O','CB','CG1','CG2'], 'LEU':['N','CA','C','O','CB','CG','CD1','CD2'], 'ILE':['N','CA','C','O','CB','CG1','CG2','CD1'], 'PRO':['N','CA','C','O','CB','CG','CD'], 'PHE':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'], 'TRP':['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'], 'MET':['N','CA','C','O','CB','CG','SD','CE'], 'GLY':['N','CA','C','O'], 'SER':['N','CA','C','O','CB','OG'], 'THR':['N','CA','C','O','CB','OG1','CG2'], 'CYS':['N','CA','C','O','CB','SG'], 'TYR':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH'], 'ASN':['N','CA','C','O','CB','CG','OD1','ND2'], 'GLN':['N','CA','C','O','CB','CG','CD','OE1','NE2'], 'ASP':['N','CA','C','O','CB','CG','OD1','OD2'], 'GLU':['N','CA','C','O','CB','CG','CD','OE1','OE2'], 'LYS':['N','CA','C','O','CB','CG','CD','CE','NZ'], 'ARG':['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'], 'HIS':['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2']}
	return x

def aminomapping():
	y={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'Q':'GLN', 'E':'GLU', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
	return y
