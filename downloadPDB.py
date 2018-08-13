import os
import urllib

path='/home/intern/Desktop/Internship/'
targetpath='/home/intern/Desktop/PDB/'

RefinedPDBCodeList = []
with open('RefinedPDBCodeList') as inputfile:
    for line in inputfile:
         RefinedPDBCodeList.append(line.rstrip('\r\n'))

for i in range(len(RefinedPDBCodeList)):
	print str(RefinedPDBCodeList[i])[:4]
	urllib.urlretrieve('http://files.rcsb.org/download/'+str(RefinedPDBCodeList[i])[:4].upper()+'.pdb', str(RefinedPDBCodeList[i]))
	os.rename(path+str(RefinedPDBCodeList[i]),targetpath+str(RefinedPDBCodeList[i]))
