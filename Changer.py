import os
import glob
import amino
import errno
from numpy import *
from scipy.optimize import *

#---------------------------------------------------------------------------------------------------------------------------------------

def fixAtom(t,L):
	x=t[0]
	y=t[1]
	z=t[2]									#function to fix the coordinates of the missing atom
										#function takes an array with direction vectors and coordinates (one each for N and CA) and an array with guess values
	F=empty((3))								#(same for both N and CA) 
										#function returns a 3D coordinate (x,y,z) where the atom is expected to be
	F[0]=(L[1]*x)-(L[0]*y)-(L[1]*L[3])+(L[0]*L[4])
	F[1]=(L[2]*y)-(L[1]*z)-(L[2]*L[4])+(L[1]*L[5])
	F[2]=(L[2]*x)-(L[0]*z)-(L[2]*L[3])+(L[0]*L[5])
	
	return F

#---------------------------------------------------------------------------------------------------------------------------------------
#initializing variables

path='/home/intern/Desktop/PDB/'
newpath='/home/intern/Desktop/UpdatedPDB/'
stdpath='/home/intern/Desktop/StandardFile/'
missing_record=[]
i=0
temp=[]
helix=[]
sheet=[]
amino_num=0
Nxco=Nyco=Nzco=0.0
CAxco=CAyco=CAzco=0.0
tempxco=tempyco=tempzco=0.0
atoms=amino.aminodict()		#a dictionary that contains a list of all the standard atoms in each amino acid

#---------------------------------------------------------------------------------------------------------------------------------------
#to extract all the ATOM records from the pdb file and storing it in another pdb file with the same name as that of the initial file in another location

for filename in glob.iglob(os.path.join(path,'*.pdb')):
	tempname=filename[-8:]
	handle=open(filename,'r')					#iterating through each of the downloaded pdb file
	lines=handle.readlines()
	newfilename=newpath+tempname

#---------------------------------------------------------------------------------------------------------------------------------------

	if not os.path.exists(os.path.dirname(newfilename)):
		try:
			os.makedirs(os.path.dirname(newfilename))	#makes sure the target directory exists
		except OSError as exc:
			if exc.errno!=errno.EEXIST:
				raise

#---------------------------------------------------------------------------------------------------------------------------------------

	with open(newfilename,'w+') as newfile:
		for eachline in lines:
			records=eachline.split()
			if records[0]=='HELIX':
				for i in range(int(eachline[21:25]),int(eachline[33:37])):	#storing the helix atoms in a list 'helix'
					helix.append(i)

			if records[0]=='SHEET':
				for i in range(int(eachline[22:26]),int(eachline[33:37])):	#storing the sheet atoms in a list 'sheet'
					sheet.append(i)
				try:
					for i in range(int(eachline[50:54]),int(eachline[65:69])):
						sheet.append(i)
				except:
					pass

			if records[0]=='ATOM':					#writing only the atom records to the new file
				chainid=eachline[21:22]
				if chainid==' ':
					chainid='A'
				if int(eachline[22:26]) in helix:
					to_write=[eachline[0:4],eachline[6:11],eachline[12:16],eachline[17:20],chainid,eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],'HELIX','\n']
					newfile.writelines('\t'.join(to_write))
				elif int(eachline[22:26]) in sheet:
					to_write=[eachline[0:4],eachline[6:11],eachline[12:16],eachline[17:20],chainid,eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],'SHEET','\n']
					newfile.writelines('\t'.join(to_write))
				else:
					to_write=[eachline[0:4],eachline[6:11],eachline[12:16],eachline[17:20],chainid,eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],'COIL','\n']
					newfile.writelines('\t'.join(to_write))
	handle.close()

#---------------------------------------------------------------------------------------------------------------------------------------
#fixing missing atoms and avoiding duplicate records if any

for filename in glob.iglob(os.path.join(newpath,'*.pdb')):	#iterating through the new pdb files which contain only ATOM records
	newfilename=newpath+'temp1.pdb'
	handle=open(filename,'r')
	newhandle=open(newfilename,'w+')
	lines=handle.readlines()
	i=0
	temp=[]							#temp is used to iterate through the dictionary
	tracker=0
	for eachline in lines:
		records=eachline.split()
		fixcoordinates1=fixcoordinates2=0.0

		if records[2]=='N' and i==len(temp):		#stores some relevant information whenever a new N atom is encountered
			amino_num=records[5]
			i=0
			temp=atoms[records[3]]			#temp contains a list of the standard atoms in the current amino acid
			aminoatomname=records[3]		#aminoatomname contains the name of the current amino acid
			tempchain=records[4]
			Nxco=float(records[6])			#stores x coordinate of N \
			Nyco=float(records[7])			#stores y coordinate of N  - for fixing missing atoms
			Nzco=float(records[8])			#stores z coordinate of N /
			helixorsheet=records[12]

		if records[2]=='CA':
			CAxco=float(records[6])			#stores x coordinate of CA \
			CAyco=float(records[7])			#stores y coordinate of CA  - for fixing missing atoms
			CAzco=float(records[8])			#stores z coordinate of CA /

		if records[5]==amino_num:			#executes if the atom is a valid member of the amino acid we've taken
			try:
				if records[2]==temp[i]:		#executes if the atom is present where it should be in the file
					to_write=[records[0],str(tracker+1),records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],helixorsheet,'\n']
					newhandle.writelines('\t'.join(to_write))
					tempxco=float(records[6])
					tempyco=float(records[7])
					tempzco=float(records[8])
					tracker+=1
					i+=1

				else:				#executes if there are atoms missing
					while records[2]!=temp[i] and records[2]!=temp[i-1]:
						#fix the missing atom

						dirfromN=[]	#to store direction vectors from N to the missing atom in the standard file
						dirfromCA=[]	#to store direction vectors from CA to the missing atom in the standard file
						paramN=[]	#parameters to be passed to the function that fixes the missing atom
						paramCA=[]	#parameters to be passed to the function that fixes the missing atom
						guess=[]	#mandatory to fix the missing atom

						stdhandle=open(stdpath+'5m5c.pdb','r')	#opens the standard file, variables starting with std
						stdlines=stdhandle.readlines()		#are used to handle the standard file

						for stdi in range(len(stdlines)):
							if stdlines[stdi][17:20].strip(' \t')!=aminoatomname:
								continue	#finding the line with the amino acid with the missing atom
							if stdlines[stdi][12:16].strip(' \t')=='N':
								dirfromN.append(float(stdlines[stdi][30:38])-float(stdlines[stdi+i][30:38]))
								dirfromN.append(float(stdlines[stdi][38:46])-float(stdlines[stdi+i][38:46]))
								dirfromN.append(float(stdlines[stdi][46:54])-float(stdlines[stdi+i][46:54]))
								dirfromCA.append(float(stdlines[stdi+1][30:38])-float(stdlines[stdi+i][30:38]))
								dirfromCA.append(float(stdlines[stdi+1][38:46])-float(stdlines[stdi+i][38:46]))
								dirfromCA.append(float(stdlines[stdi+1][46:54])-float(stdlines[stdi+i][46:54]))
								break	#after storing the direction vectors, exiting the loop
						stdhandle.close()

						for eachvalue in dirfromN:
							paramN.append(eachvalue)#\
						paramN.append(Nxco)		 #\
						paramN.append(Nyco)		 # - populating the parameter list for N 
						paramN.append(Nzco)		 #/

						for eachvalue in dirfromCA:
							paramCA.append(eachvalue)#\
						paramCA.append(CAxco)		  #\
						paramCA.append(CAyco)		  # - populating the parameter list for CA
						paramCA.append(CAzco)		  #/

						guess=[(tempxco+float(records[6]))/2,(tempyco+float(records[7]))/2,(tempzco+float(records[8]))/2]	#calculating a guess value based on the coordinates of the atoms above and below the missing atom

						fixcoordinates1=fsolve(fixAtom,guess,paramN)	#x,y,z values from N
						fixcoordinates2=fsolve(fixAtom,guess,paramCA)	#x,y,z values from CA

						finalx=(fixcoordinates1[0]+fixcoordinates2[0])/2 #\
						finaly=(fixcoordinates1[1]+fixcoordinates2[1])/2 # - x,y,z coordinates of the missing atom
						finalz=(fixcoordinates1[2]+fixcoordinates2[2])/2 #/

						missing_record=[records[0],str(tracker+1),temp[i],aminoatomname,records[4],records[5],str(finalx)[:6],str(finaly)[:6],str(finalz)[:6],'1.00','99.99',temp[i][0],helixorsheet,'\n']
						newhandle.writelines('\t'.join(missing_record))
						tempxco=float(str(finalx)[:6])
						tempyco=float(str(finaly)[:6])
						tempzco=float(str(finalz)[:6])
						tracker+=1
						i+=1

					if records[2]!=temp[i-1]:	#skipping duplicate record, if present
						to_write=[records[0],str(tracker+1),records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],helixorsheet,'\n']
						newhandle.writelines('\t'.join(to_write))
						tempxco=float(records[6])
						tempyco=float(records[7])
						tempzco=float(records[8])
						tracker+=1
						i+=1

			except:				#to write atoms other than the standard atoms in the amino acid, if any
				to_write=[records[0],str(tracker+1),records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],helixorsheet,'\n']
				newhandle.writelines('\t'.join(to_write))
				tempxco=float(records[6])
				tempyco=float(records[7])
				tempzco=float(records[8])
				tracker+=1

		else:
			
			while i<len(temp):
				#fix the missing atom

				dirfromN=[]		#to store direction vectors from N to the missing atom in the standard file
				dirfromCA=[]		#to store direction vectors from CA to the missing atom in the standard file
				paramN=[]		#parameters to be passed to the function that fixes the missing atom
				paramCA=[]		#parameters to be passed to the function that fixes the missing atom
				guess=[]		#mandatory to fix the missing atom

				stdhandle=open(stdpath+'5m5c.pdb','r')	#opens the standard file, variables starting with std
				stdlines=stdhandle.readlines()		#are used to handle the standard file

				for stdi in range(len(stdlines)):
					if stdlines[stdi][17:20].strip(' \t')!=aminoatomname:
						continue		#finding the line with the amino acid with the missing atom
					if stdlines[stdi][12:16].strip(' \t')=='N':
						dirfromN.append(float(stdlines[stdi][30:38])-float(stdlines[stdi+i][30:38]))
						dirfromN.append(float(stdlines[stdi][38:46])-float(stdlines[stdi+i][38:46]))
						dirfromN.append(float(stdlines[stdi][46:54])-float(stdlines[stdi+i][46:54]))
						dirfromCA.append(float(stdlines[stdi+1][30:38])-float(stdlines[stdi+i][30:38]))
						dirfromCA.append(float(stdlines[stdi+1][38:46])-float(stdlines[stdi+i][38:46]))
						dirfromCA.append(float(stdlines[stdi+1][46:54])-float(stdlines[stdi+i][46:54]))
						break		#after storing the direction vectors, exiting the loop
				stdhandle.close()

				for eachvalue in dirfromN:
					paramN.append(eachvalue)#\
				paramN.append(Nxco)		 #\
				paramN.append(Nyco)		 # - populating the parameter list for N
				paramN.append(Nzco)		 #/

				for eachvalue in dirfromCA:
					paramCA.append(eachvalue)#\
				paramCA.append(CAxco)		  #\
				paramCA.append(CAyco)		  # - populating the parameter list for CA
				paramCA.append(CAzco)		  #/

					#calculating a guess value based on the coordinates of the atoms above and below the missing atom
				guess=[(tempxco+float(records[6]))/2,(tempyco+float(records[7]))/2,(tempzco+float(records[8]))/2]

				fixcoordinates1=fsolve(fixAtom,guess,paramN)	#x,y,z values from N
				fixcoordinates2=fsolve(fixAtom,guess,paramCA)	#x,y,z values from CA

				finalx=(fixcoordinates1[0]+fixcoordinates2[0])/2	#\
				finaly=(fixcoordinates1[1]+fixcoordinates2[1])/2	# - x,y,z coordinates of the missing atom
				finalz=(fixcoordinates1[2]+fixcoordinates2[2])/2	#/

				missing_record=[records[0],str(tracker+1),temp[i],aminoatomname,tempchain,amino_num,str(finalx)[:6],str(finaly)[:6],str(finalz)[:6],'1.00','99.99',temp[i][0],helixorsheet,'\n']
				newhandle.writelines('\t'.join(missing_record))
				tempxco=float(str(finalx)[:6])
				tempyco=float(str(finaly)[:6])
				tempzco=float(str(finalz)[:6])
				i+=1
				tracker+=1

			

			to_write=[records[0],str(tracker+1),records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],helixorsheet,'\n']
			newhandle.writelines('\t'.join(to_write))
			tempxco=float(records[6])
			tempyco=float(records[7])
			tempzco=float(records[8])
			tracker+=1

			amino_num=records[5]
			i=1
			temp=atoms[records[3]]
			aminoatomname=records[3]
			tempchain=records[4]
			Nxco=float(records[6])
			Nyco=float(records[7])
			Nzco=float(records[8])
			helixorsheet=records[12]

	handle.close()
	newhandle.close()
	os.rename(newfilename,newpath+filename[-8:])			#replacing the input file with a file that has fixed the missing atoms and duplicate records, if any

#---------------------------------------------------------------------------------------------------------------------------------------
