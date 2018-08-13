import os
import glob
import amino
import errno
from numpy import *
from scipy.optimize import *
from prody import *
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
path='/home/jisna_mary/Desktop/'
libpath='/home/jisna_mary/Desktop/Internship_new/'
lib2path='/home/jisna_mary/Desktop/Internship_new/Library2/'
targetpath='/home/jisna_mary/Desktop/Internship_new/Run/'
stdpath='/home/jisna_mary/Desktop/Internship_new/StandardFile/'
pymolpath='/home/jisna_mary/Desktop/Internship_new/pymol/'

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Function to change files to pdb format

def writePDBFormat(frompath,frompdb):
	tracker=1
	aminonum=0
	temppdb=open(frompath,'w+')	#frompath=path+'temp.pdb'
	pdbfile=open(frompdb,'r')	#frompdb=path+'WXYZ.pdb'
	pdblines=pdbfile.readlines()
	pdbfile.close()	
	for eachline in pdblines:
		records=eachline.split()
		if records[2].strip(' \t')=='N':
			aminonum+=1
		records[0]=records[0].ljust(6)					#atom
		records[1]=str(tracker).rjust(5)				#atomnum
		records[2]=records[2].center(4)					#atomname
		records[3]=records[3].ljust(3)					#resname
		records[4]=chain.rjust(1)					#Astring
		records[5]=str(aminonum).rjust(4)				#resnum
		records[6]=str('%8.3f' % (float(records[6]))).rjust(8)		#x
		records[7]=str('%8.3f' % (float(records[7]))).rjust(8)		#y
		records[8]=str('%8.3f' % (float(records[8]))).rjust(8)		#z
		records[9]=str('%6.2f' % (float(records[9]))).rjust(6)		#occ
		records[10]=str('%6.2f' % (float(records[10]))).ljust(6)	#temp
		records[11]=records[11].rjust(12)				#elemname    
		temppdb.write('%s%s %s %s %s%s    %s%s%s%s%s%s\n' % (records[0],records[1],records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11]))
		tracker+=1

	if aminonum<input_fasta_len:
		print '\nThere is(are) ',input_fasta_len-aminonum,' missing residue(s) in the output\n'

	temppdb.close()
	os.rename(frompath,pymolpath+'WXYZ.pdb')

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Initialization block

mapping=amino.aminomapping()
atoms=amino.aminodict()
subject=[]
query=[]
temp=[]
clusters=[]
helixorsheet=''
first_match=-1
model_index=0
smallestyindex=0
amino_num=0
tracker=0
slno=1
flag=0
pos=0
j=k=l=0
Nxco=Nyco=Nzco=0.0
CAxco=CAyco=CAzco=0.0
tempxco=tempyco=tempzco=0.0
xco=yco=zco=0.0

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

input_fasta=raw_input('Enter the FASTA sequence:\n')			#user inputs the FASTA sequence
input_fasta_len=len(input_fasta)
mini=input_fasta_len
f = open('1k8k.fasta.txt', 'r+')
header=f.readline()
f.seek(0)
f.truncate()
f.close()
f = open('1k8k.fasta.txt', 'w')
f.write(header)
for i in range(0,input_fasta_len+1,80):
	f.writelines(input_fasta[i:i+80])				#writes the input sequence according to the FASTA file format (80 amino acids in each line)
	f.write('\n')
f.close()

print '\nLength of the input FASTA sequence is\t',input_fasta_len,'\n'

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for i in range(input_fasta_len):				#converting one letter code to three letter code
	query.append(mapping[input_fasta[i]])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#BLASTing

print '\nProcessing\n'

f_record = next(SeqIO.parse('1k8k.fasta.txt', 'fasta'))
my_query = SeqIO.read("1k8k.fasta.txt", format="fasta") 
result_handle = NCBIWWW.qblast("blastp", "pdb", f_record.format('fasta'))	# if BLASTs for homo sapiens -> ,entrez_query='HOMO SAPIENS'
blast_result = open("my_blast.xml", "w")
blast_result.write(result_handle.read())
blast_result.close()
result_handle.close()

use_handle=open("my_blast.xml")				#opens the XML file which contains the blast results
blast_records=NCBIXML.parse(use_handle)
blast_record=blast_records.next()

try:
	alignment= blast_record.alignments[0]
except IndexError:
	print '\nNo match after Blasting.\n'
	exit()
title= alignment.title

chain=title.split("|")[4][0]
query_from=blast_record.alignments[0].hsps[0].query_start
query_to=blast_record.alignments[0].hsps[0].query_end
hit_from=blast_record.alignments[0].hsps[0].sbjct_start
hit_to=blast_record.alignments[0].hsps[0].sbjct_end

match=blast_record.alignments[0].hsps[0].match
subject_len=len(match)
for i in range(subject_len):
	if match[i]=='+' or match[i]==' ':
		subject.append('_')
	else:
		subject.append(mapping[match[i]])
print blast_record.alignments[0].hsps[0].sbjct
print 'Subject Length= ',subject_len
print blast_record.alignments[0].hsps[0].query

protein_id= title.split("|")[3]				#Extracts the pdb file's name of the first match
print '\nProtein Code\t',protein_id.lower()
results=fetchPDB(protein_id,compressed=False)

os.rename(libpath+protein_id.lower()+'.pdb',targetpath+protein_id.lower()+'.pdb')	#Moves the file to another location
print '\n'+protein_id.lower()+'.pdb moved to location '+targetpath+'\n'

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

intermediate_file=open(path+protein_id.lower()+'.pdb','w+')
with open(targetpath+protein_id.lower()+'.pdb') as pdbfile:
	lines=pdbfile.readlines()
	for eachline in lines:
		records=eachline.split()
		if records[0]=='ATOM':
			if eachline[12:16].strip(' \t')=='N':
				intermediate_file.write(eachline)
			if eachline[12:16].strip(' \t')=='CA':
				intermediate_file.write(eachline)
		elif records[0]=='SEQRES':
			records.append('\n')
			intermediate_file.writelines(' '.join(records))
		elif records[0]=='TER':
			intermediate_file.write(eachline)
			slno+=1
		elif records[0]=='REMARK':		#REMARK 465 records contain the missing residues
			if records[1]=='465' and len(records)==5:
				records.append('\n')
				intermediate_file.writelines('\t'.join(records))
		else:
			continue
intermediate_file.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fixing missing atoms and eliminating duplicate records in the downloaded pdb file

intermediate_file=open(targetpath+protein_id.lower()+'.pdb','w+')
with open(path+protein_id.lower()+'.pdb') as pdbfile:
	i=0
	slno=1
	temp=[]
	lines=pdbfile.readlines()
	for eachline in lines:
		records=eachline.split()

		if records[0]=='ATOM':
			if records[2]=='N' and i==len(temp):		#stores some relevant information whenever a new N atom is encountered
				amino_num=int(eachline[22:26])
				i=0
				temp=atoms[eachline[17:20]]		#temp contains a list of the standard atoms in the current amino acid
				aminoatomname=eachline[17:20]		#aminoatomname contains the name of the current amino acid
				tempchain=eachline[21:22]
				Nxco=float(eachline[30:38])			#stores x coordinate of N \
				Nyco=float(eachline[38:46])			#stores y coordinate of N  - for fixing missing atoms
				Nzco=float(eachline[46:54])			#stores z coordinate of N /
				try:
					helixorsheet=records[12]
				except:
					helixorsheet=''

			if records[2]=='CA':
				CAxco=float(eachline[30:38])			#stores x coordinate of CA \
				CAyco=float(eachline[38:46])			#stores y coordinate of CA  - for fixing missing atoms
				CAzco=float(eachline[46:54])			#stores z coordinate of CA /

			if int(eachline[22:26])==amino_num:		#executes if the atom is a valid member of the amino acid we've taken
				try:
					if eachline[12:16].strip(' \t')==temp[i]:	#executes if the atom is present where it should be
						to_write=[eachline[0:4],str(tracker+1),eachline[12:16],eachline[17:20],eachline[21:22],eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],helixorsheet,'\n']
						intermediate_file.writelines('\t'.join(to_write))
						tempxco=float(eachline[30:38])
						tempyco=float(eachline[38:46])
						tempzco=float(eachline[46:54])
						tracker+=1
						i+=1

					else:							#executes if there are atoms missing
						while eachline[12:16].strip(' \t')!=temp[i] and eachline[12:16].strip(' \t')!=temp[i-1]:
							dirfromN=[]			#to store direction vectors from N to the missing atom in the standard file
							dirfromCA=[]			#to store direction vectors from CA to the missing atom in the standard file
							paramN=[]			#parameters to be passed to the function that fixes the missing atom
							paramCA=[]			#parameters to be passed to the function that fixes the missing atom
							guess=[]			#mandatory to fix the missing atom

							stdhandle=open(stdpath+'5d8v.pdb','r')	#opens the standard file, variables starting with std are used to handle the standard file
							stdlines=stdhandle.readlines()

							for stdi in range(len(stdlines)):
								if stdlines[stdi][17:20].strip(' \t')!=aminoatomname:
									continue     			#finding the line with the amino acid with the missing atom
								if stdlines[stdi][12:16].strip(' \t')=='N':
									dirfromN.append(float(stdlines[stdi][30:38])-float(stdlines[stdi+i][30:38]))
									dirfromN.append(float(stdlines[stdi][38:46])-float(stdlines[stdi+i][38:46]))
									dirfromN.append(float(stdlines[stdi][46:54])-float(stdlines[stdi+i][46:54]))
									dirfromCA.append(float(stdlines[stdi+1][30:38])-float(stdlines[stdi+i][30:38]))
									dirfromCA.append(float(stdlines[stdi+1][38:46])-float(stdlines[stdi+i][38:46]))
									dirfromCA.append(float(stdlines[stdi+1][46:54])-float(stdlines[stdi+i][46:54]))
									break			#after storing the direction vectors, exiting the loop
							stdhandle.close()

							for eachvalue in dirfromN:
								paramN.append(eachvalue)	#\
							paramN.append(Nxco)			 #\
							paramN.append(Nyco)		 	# - populating the parameter list for N 
							paramN.append(Nzco)		 	#/

							for eachvalue in dirfromCA:
								paramCA.append(eachvalue)	#\
							paramCA.append(CAxco)		 	 #\
							paramCA.append(CAyco)		 	 # - populating the parameter list for CA
							paramCA.append(CAzco)		 	 #/

							guess=[(tempxco+float(records[6]))/2,(tempyco+float(records[7]))/2,(tempzco+float(records[8]))/2]	#calculating a guess value based on the coordinates of the atoms above and below the missing atom

							fixcoordinates1=fsolve(fixAtom,guess,paramN)	#x,y,z values from N
							fixcoordinates2=fsolve(fixAtom,guess,paramCA)	#x,y,z values from CA

							finalx=(fixcoordinates1[0]+fixcoordinates2[0])/2 	#\
							finaly=(fixcoordinates1[1]+fixcoordinates2[1])/2 	#-x,y,z coordinates of the missing atom
							finalz=(fixcoordinates1[2]+fixcoordinates2[2])/2 	#/

							missing_record=[records[0],str(tracker+1),temp[i],aminoatomname,records[4],records[5],str(finalx)[:6],str(finaly)[:6],str(finalz)[:6],'1.00','99.99',temp[i][0],helixorsheet,'\n']
							intermediate_file.writelines('\t'.join(missing_record))
							tempxco=float(str(finalx)[:6])
							tempyco=float(str(finaly)[:6])
							tempzco=float(str(finalz)[:6])
							tracker+=1
							i+=1

						if eachline[12:16].strip(' \t')!=temp[i-1]:	#skipping duplicate records, if present
							to_write=[eachline[0:4],str(tracker+1),eachline[12:16],eachline[17:20],eachline[21:22],eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],helixorsheet,'\n']
							intermediate_file.writelines('\t'.join(to_write))
							tempxco=float(eachline[30:38])
							tempyco=float(eachline[38:46])
							tempzco=float(eachline[46:54])
							tracker+=1
							i+=1

				except:				#to write atoms other than the standard atoms in the amino acid, if any
					if eachline[12:16].strip(' \t')!=temp[i-1]:
						to_write=[eachline[0:4],str(tracker+1),eachline[12:16],eachline[17:20],eachline[21:22],eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],helixorsheet,'\n']
						intermediate_file.writelines('\t'.join(to_write))
						tempxco=float(eachline[30:38])
						tempyco=float(eachline[38:46])
						tempzco=float(eachline[46:54])
						tracker+=1

			else:
				while i<len(temp):
					dirfromN=[]		#to store direction vectors from N to the missing atom in the standard file
					dirfromCA=[]		#to store direction vectors from CA to the missing atom in the standard file
					paramN=[]		#parameters to be passed to the function that fixes the missing atom
					paramCA=[]		#parameters to be passed to the function that fixes the missing atom
					guess=[]		#mandatory to fix the missing atom

					stdhandle=open(stdpath+'5d8v.pdb','r')	#opens the standard file, variables starting with std are used to handle the standard file
					stdlines=stdhandle.readlines()

					for stdi in range(len(stdlines)):
						if stdlines[stdi].split()[0]=='ATOM':
							if stdlines[stdi][17:20].strip(' \t')!=aminoatomname:
								continue			#finding the line with the amino acid with the missing atom
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
						paramN.append(eachvalue)	#\
					paramN.append(Nxco)			 #\
					paramN.append(Nyco)			 # - populating the parameter list for N
					paramN.append(Nzco)			 #/

					for eachvalue in dirfromCA:
						paramCA.append(eachvalue)	#\
					paramCA.append(CAxco)		 	 #\
					paramCA.append(CAyco)		 	 # - populating the parameter list for CA
					paramCA.append(CAzco)		 	 #/

					#calculating a guess value based on the coordinates of the atoms above and below the missing atom
					guess=[(tempxco+float(records[6]))/2,(tempyco+float(records[7]))/2,(tempzco+float(records[8]))/2]

					fixcoordinates1=fsolve(fixAtom,guess,paramN)	#x,y,z values from N
					fixcoordinates2=fsolve(fixAtom,guess,paramCA)	#x,y,z values from CA

					finalx=(fixcoordinates1[0]+fixcoordinates2[0])/2   	#\
					finaly=(fixcoordinates1[1]+fixcoordinates2[1])/2   	# - x,y,z coordinates of the missing atom
					finalz=(fixcoordinates1[2]+fixcoordinates2[2])/2   	#/

					missing_record=[records[0],str(tracker+1),temp[i],aminoatomname,tempchain,str(amino_num),str(finalx)[:6],str(finaly)[:6],str(finalz)[:6],'1.00','99.99',temp[i][0],helixorsheet,'\n']
					intermediate_file.writelines('\t'.join(missing_record))
					tempxco=float(str(finalx)[:6])
					tempyco=float(str(finaly)[:6])
					tempzco=float(str(finalz)[:6])
					i+=1
					tracker+=1				

				to_write=[eachline[0:4],str(tracker+1),eachline[12:16],eachline[17:20],eachline[21:22],eachline[22:27],eachline[30:38],eachline[38:46],eachline[46:54],eachline[54:60],eachline[60:66],eachline[72:76],eachline[76:78],helixorsheet,'\n']
				intermediate_file.writelines('\t'.join(to_write))
				tempxco=float(eachline[30:38])
				tempyco=float(eachline[38:46])
				tempzco=float(eachline[46:54])
				tracker+=1

				amino_num=int(eachline[22:26])
				i=1
				temp=atoms[eachline[17:20]]
				aminoatomname=eachline[17:20]
				tempchain=eachline[21:22]
				Nxco=float(eachline[30:38])
				Nyco=float(eachline[38:46])
				Nzco=float(eachline[46:54])
				try:
					helixorsheet=records[12]
				except:
					helixorsheet=''

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		elif records[0]=='SEQRES':
			records.append('\n')
			intermediate_file.writelines(' '.join(records))
		elif records[0]=='TER':
			intermediate_file.write(eachline)
			slno+=1
		elif records[0]=='REMARK':		#REMARK 465 records contain the missing residues
			if records[1]=='465' and len(records)==5:
				records.append('\n')
				intermediate_file.writelines('\t'.join(records))
		else:
			continue
intermediate_file.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filling the gaps in the match

j=k=0
first_match=query_from-1
newpdb=open(path+'WXYZ.pdb','w+')
thandle=open(path+'temp_model.txt','w+')
thandle.close()
usepdb=open(targetpath+protein_id.lower()+'.pdb','r')
usepdblines=usepdb.readlines()
tracker=0
aminonum=1

for i in range(len(usepdblines)):
	records=usepdblines[i].split()
	if records[0]=='REMARK' and int(usepdblines[i][17:20])<=0 and records[3]==chain:	#adjusts the values of hit_from and hit_to for every residue shown with a non-positive value under SSEQI
		hit_from-=1
		hit_to-=1
	if records[0]!='REMARK':
		break

m=hit_from
n=0
while m<=hit_to:
	n=0
	while True:
		records=usepdblines[n].split()
		if records[0]=='REMARK' and int(usepdblines[n][17:20])==m and records[3]==chain:	#alters the subject list if one of our residues is present under the REMARK 465 section
			subject[int(usepdblines[n][17:20])-hit_from]='_'
		n+=1
		if records[0]!='REMARK':
			break
	m+=1

m=0
while True:
	if subject[m]=='_':	#if there are underscores at the start of the subject list, first_match, k, hit_from and hit_to are adjusted
		first_match+=1
		hit_from+=1
		hit_to+=1
		k+=1		#k is used to iterate through the subject list, the value of k has to be properly maintained
	if subject[m]!='_':
		break
	m+=1

m=hit_from
j=first_match
l=first_match
print 'First match ',first_match

if first_match==0:			#executes if there are no gaps to the left of the first match
	j=first_match
	k=first_match
	l=first_match
	while k<subject_len:				#starts execution after the gaps to the left are filled, iterates till the end of the query sequence
		if blast_record.alignments[0].hsps[0].query[l-first_match]=='-':	#skips through the hyphens in the query sequence of the BLAST result
			l+=1
			k+=1
			m+=1
			continue
		
		try:
			query[j]
		except:
			newpdb.close()
			writePDBFormat(path+'temp.pdb',path+'WXYZ.pdb')
			exit()

		if query[j]==subject[k]:
			for eachindex in range(len(usepdblines)):
				records=usepdblines[eachindex].split()
				if records[0]!='ATOM':
					continue
				if records[4]==chain and m==int(records[5]) and records[3]==query[j]:		#copies from the pdb file if the query matches the blast result
					a=eachindex
					try:
						while True:
							records=usepdblines[a].split()
							if int(records[5])!=m:
								break
							if records[2]=='CA':
								xco=float(records[6])
								yco=float(records[7])
								zco=float(records[8])
							to_write=[records[0],str(tracker+1),records[2],records[3],records[4],str(aminonum),records[6],records[7],records[8],records[9],records[10],records[11],'\n']
							newpdb.write('\t'.join(to_write))
							tracker+=1
							a+=1
						break

					except:
						print query[j],
						aminonum+=1
						j+=1
						while j<input_fasta_len:
							libfile=open(lib2path+query[j-1]+'_'+query[j]+'.txt','r')
							libfilelines=libfile.readlines()
							clusterline=libfilelines[0]
							smallest=99999999
							for i in range(len(clusterline.split())):
								if i!=0:
									clusters.append(float(clusterline.split()[i]))
							for each in clusters:
								if abs(each-xco)<abs(smallest-xco):
									smallest=float(each)			#identifies the cluster closest to the xco of the last CA

							tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file
							tempcoords.close()

							for eachindex in range(len(libfilelines)):			#From the library file, y co-ordinates of all the models in the closest cluster are copied to a temp file so that we can
								records=libfilelines[eachindex].split()	#select the model that has the closest value for the y co-ordinate when compared to the value stored in yco
								a=eachindex
								tempcoords=open(path+'tempcoords.txt','a')
								if records[0]=='MODEL' and float(records[2][:7])==smallest:
									while libfilelines[a].split()[0]!='ENDMODEL':
										if libfilelines[a].split()[0]=='MODEL':
											tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
										a+=1
								tempcoords.close()

								smallesty=999999999
								smallestyindex=-1
								tempcoords=open(path+'tempcoords.txt','r')
								templines=tempcoords.readlines()
								for eachrecord in templines:				#Finds the model number of the model with the closest value for y co-ordinate
									splitrec=eachrecord.split()
									if abs(float(splitrec[1])-yco)<abs(smallesty):
										smallesty=float(splitrec[1])-yco
										smallestyindex=splitrec[0]
								tempcoords.close()

							for eachindex in range(len(libfilelines)):		#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
								records=libfilelines[eachindex].split()		#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
								if records[0]=='MODEL':				#closest to the value in yco
									if float(records[2][:7])==smallest:
										if records[1]==smallestyindex:
											a=eachindex
											while True:
												libfilerecords=libfilelines[a].split()
												if libfilerecords[0]=='MODEL':
													pass
												elif libfilerecords[0]=='ENDMODEL':
													break
												else:
													if libfilerecords[1]=='N':
														a+=1
														while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file
															a+=1
													libfilerecords=libfilelines[a].split()
													to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
													if libfilerecords[1]=='CA':
														xco+=float(libfilerecords[4])
														yco+=float(libfilerecords[5])
														zco+=float(libfilerecords[6])
													newpdb.write('\t'.join(to_write))
													tracker+=1
												a+=1
											break
							libfile.close()
							aminonum+=1
							print '_',
							j+=1
							clusters=[]
						newpdb.close()
						writePDBFormat(path+'temp.pdb',path+'WXYZ.pdb')
						exit()
			aminonum+=1
			print query[j],
			j+=1
			k+=1
			l+=1
			m+=1
		else:
			libfile=open(lib2path+query[j-1]+'_'+query[j]+'.txt','r')
			libfilelines=libfile.readlines()
			clusterline=libfilelines[0]
			smallest=99999999
			for i in range(len(clusterline.split())):
				if i!=0:
					clusters.append(float(clusterline.split()[i]))
			for each in clusters:
				if abs(each-xco)<abs(smallest-xco):
					smallest=float(each)			#identifies the cluster closest to the xco of the last CA

			tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file
			tempcoords.close()

			for eachindex in range(len(libfilelines)):			#From the library file, y co-ordinates of all the models in the closest cluster are copied to a temp file so that we can
				records=libfilelines[eachindex].split()			#select the model that has the closest value for the y co-ordinate when compared to the value stored in yco
				a=eachindex
				tempcoords=open(path+'tempcoords.txt','a')
				if records[0]=='MODEL' and float(records[2][:7])==smallest:
					while libfilelines[a].split()[0]!='ENDMODEL':
						if libfilelines[a].split()[0]=='MODEL':
							tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
						a+=1
				tempcoords.close()

				smallesty=999999999
				smallestyindex=-1
				tempcoords=open(path+'tempcoords.txt','r')
				templines=tempcoords.readlines()
				for eachrecord in templines:				#Finds the model number of the model with the closest value for y co-ordinate
					splitrec=eachrecord.split()
					if abs(float(splitrec[1])-yco)<abs(smallesty):
						smallesty=float(splitrec[1])-yco
						smallestyindex=splitrec[0]
				tempcoords.close()

			for eachindex in range(len(libfilelines)):		#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
				records=libfilelines[eachindex].split()		#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
				if records[0]=='MODEL':				#closest to the value in yco
					if float(records[2][:7])==smallest:
						if records[1]==smallestyindex:
							a=eachindex
							while True:
								libfilerecords=libfilelines[a].split()
								if libfilerecords[0]=='MODEL':
									pass
								elif libfilerecords[0]=='ENDMODEL':
									break
								else:
									if libfilerecords[1]=='N':
										a+=1
										while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file as we have already written it
											a+=1
									libfilerecords=libfilelines[a].split()
									to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
									if libfilerecords[1]=='CA':
										xco+=float(libfilerecords[4])
										yco+=float(libfilerecords[5])
										zco+=float(libfilerecords[6])
									newpdb.write('\t'.join(to_write))
									tracker+=1
								a+=1
							break
			aminonum+=1
			libfile.close()
			clusters=[]
			print '_',
			if blast_record.alignments[0].hsps[0].sbjct[k-first_match]=='-':	#skips through the hyphens in the subject sequence of the BLAST result
				j+=1
				k+=1
				l+=1
			else:
				j+=1
				k+=1
				l+=1
				m+=1

	while j<input_fasta_len:
		libfile=open(lib2path+query[j-1]+'_'+query[j]+'.txt','r')
		libfilelines=libfile.readlines()
		clusterline=libfilelines[0]
		smallest=99999999
		for i in range(len(clusterline.split())):
			if i!=0:
				clusters.append(float(clusterline.split()[i]))
		for each in clusters:
			if abs(each-xco)<abs(smallest-xco):
				smallest=float(each)			#identifies the cluster closest to the xco of the last CA

		tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file
		tempcoords.close()

		for eachindex in range(len(libfilelines)):			#From the library file, y co-ordinates of all the models in the closest cluster are copied to a temp file so that we can
			records=libfilelines[eachindex].split()			#select the model that has the closest value for the y co-ordinate when compared to the value stored in yco
			a=eachindex
			tempcoords=open(path+'tempcoords.txt','a')
			if records[0]=='MODEL' and float(records[2][:7])==smallest:
				while libfilelines[a].split()[0]!='ENDMODEL':
					if libfilelines[a].split()[0]=='MODEL':
						tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
					a+=1
			tempcoords.close()

			smallesty=999999999
			smallestyindex=-1
			tempcoords=open(path+'tempcoords.txt','r')
			templines=tempcoords.readlines()
			for eachrecord in templines:				#Finds the model number of the model with the closest value for y co-ordinate
				splitrec=eachrecord.split()
				if abs(float(splitrec[1])-yco)<abs(smallesty):
					smallesty=float(splitrec[1])-yco
					smallestyindex=splitrec[0]
			tempcoords.close()

		for eachindex in range(len(libfilelines)):		#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
			records=libfilelines[eachindex].split()		#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
			if records[0]=='MODEL':				#closest to the value in yco
				if float(records[2][:7])==smallest:
					if records[1]==smallestyindex:
						a=eachindex
						while True:
							libfilerecords=libfilelines[a].split()
							if libfilerecords[0]=='MODEL':
								pass
							elif libfilerecords[0]=='ENDMODEL':
								break
							else:
								if libfilerecords[1]=='N':
									a+=1
									while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file as we have already written it
										a+=1
								libfilerecords=libfilelines[a].split()
								to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
								if libfilerecords[1]=='CA':
									xco+=float(libfilerecords[4])
									yco+=float(libfilerecords[5])
									zco+=float(libfilerecords[6])
								newpdb.write('\t'.join(to_write))
								tracker+=1
							a+=1
						break
		aminonum+=1
		libfile.close()
		clusters=[]
		print '_',
		j+=1
		m+=1

else:					#executes if there are gaps to the left of the first match
	m=hit_from
	for eachline in usepdblines:
		records=eachline.split()
		#print records,'\n'
		#continue
		if records[0]!='ATOM':
			continue
		if records[4]==chain and m==int(records[5]):
			if records[2]=='CA':
				xco=float(records[6])	#\
				yco=float(records[7])	# - Stores the co-ordinate values of the first match because it is used to fix the gaps to the left of the first match
				zco=float(records[8])	#/

	b=first_match
	j=b-1
	while j>=0:
		thandle=open(path+'temp_model.txt','a')
		libfile=open(lib2path+query[b]+'_'+query[j]+'.txt','r')
		libfilelines=libfile.readlines()
		clusterline=libfilelines[0]
		smallest=9999999
		for i in range(len(clusterline.split())):
			if i!=0:
				clusters.append(float(clusterline.split()[i]))
		for each in clusters:
			if abs(each-xco)<abs(smallest-xco):		#identifies the cluster closest to the xco of the last CA
				smallest=float(each)

		tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file	
		tempcoords.close()

		for eachindex in range(len(libfilelines)):			#From the library file, y co-ordinates of all the models in the	closest cluster are copied to a temp file so that we can select
			records=libfilelines[eachindex].split()			#the model that has the closest value for the y co-ordinate when compared to the value stored in yco
			a=eachindex
			tempcoords=open(path+'tempcoords.txt','a')
			if records[0]=='MODEL' and float(records[2][:7])==smallest:
				while libfilelines[a].split()[0]!='ENDMODEL':
					if libfilelines[a].split()[0]=='MODEL':
						tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
					a+=1
			tempcoords.close()

			smallesty=999999999
			smallestyindex=-1
			tempcoords=open(path+'tempcoords.txt','r')
			templines=tempcoords.readlines()
			for eachrecord in templines:					#Finds the model number of the model with the closest value for y co-ordinate
				splitrec=eachrecord.split()
				if abs(float(splitrec[1])-yco)<abs(smallesty):
					smallesty=float(splitrec[1])-yco
					smallestyindex=splitrec[0]
			tempcoords.close()

		for eachindex in range(len(libfilelines)):			#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
			records=libfilelines[eachindex].split()			#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
			if records[0]=='MODEL':					#closest to the value in yco
				if float(records[2][:7])==smallest:
					if records[1]==smallestyindex:
						a=eachindex
						while True:
							libfilerecords=libfilelines[a].split()
							if libfilerecords[0]=='MODEL':
								thandle.write('MODEL\n')
							elif libfilerecords[0]=='ENDMODEL':
								break
							else:
								if libfilerecords[1]=='N':
									a+=1
									while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file as we have already written it
										a+=1
								libfilerecords=libfilelines[a].split()
								to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
								if libfilerecords[1]=='CA':
									xco+=float(libfilerecords[4])
									yco+=float(libfilerecords[5])
									zco+=float(libfilerecords[6])
								thandle.write('\t'.join(to_write))	#The model is written to a temp file because the order of the amino acids has to be reversed
								tracker+=1
							a+=1
						thandle.write('ENDMODEL\n')
						break
		aminonum+=1
		thandle.close()
		libfile.close()
		clusters=[]
		print '_',
		b-=1
		j=b-1

	j=first_match
	thandle=open(path+'temp_model.txt','r')
	tlines=thandle.readlines()
	while model_index!=1:
		model_index=0
		for eachtindex in range(len(tlines)):				#From the temp file, the models are written to the new pdb file in the correct order (the correct order is the reverse of the order 
			trecords=tlines[eachtindex].split()			#of the models in the temp file)
			if j>0:
				if trecords[0]=='MODEL':
					model_index+=1
				if model_index==j:
					x=eachtindex
					while tlines[x].split()[0]!='ENDMODEL':
						if tlines[x].split()[0]!='MODEL':
							newpdb.writelines(tlines[x])
						x+=1
					j-=1
			else:
				break
	thandle.close()

	j=first_match
	m=hit_from
	while k<subject_len:				#starts execution after the gaps to the left are filled, iterates till the end of the query sequence
		if blast_record.alignments[0].hsps[0].query[l-first_match]=='-':	#skip through the hyphens in the query sequence of the BLAST result
			l+=1
			k+=1
			m+=1
			continue
		if query[j]==subject[k]:
			for eachindex in range(len(usepdblines)):
				records=usepdblines[eachindex].split()
				if records[0]!='ATOM':
					continue
				if records[4]==chain and m==int(records[5]) and records[3]==query[j]:		#copies from the pdb file if the query matches the blast result
					a=eachindex
					try:
						while True:
							records=usepdblines[a].split()
							if int(records[5])!=m:
								break
							if records[2]=='CA':
								xco=float(records[6])
								yco=float(records[7])
								zco=float(records[8])
							to_write=[records[0],str(tracker+1),records[2],records[3],records[4],str(aminonum),records[6],records[7],records[8],records[9],records[10],records[11],'\n']
							newpdb.write('\t'.join(to_write))
							tracker+=1
							a+=1
						break

					except:
						print query[j],
						aminonum+=1
						j+=1
						while j<input_fasta_len:
							libfile=open(lib2path+query[j-1]+'_'+query[j]+'.txt','r')
							libfilelines=libfile.readlines()
							clusterline=libfilelines[0]
							smallest=99999999
							for i in range(len(clusterline.split())):
								if i!=0:
									clusters.append(float(clusterline.split()[i]))
							for each in clusters:
								if abs(each-xco)<abs(smallest-xco):
									smallest=float(each)			#identifies the cluster closest to the xco of the last CA

							tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file
							tempcoords.close()

							for eachindex in range(len(libfilelines)):#From the library file, y co-ordinates of all the models in the closest cluster are copied to a temp file so that
								records=libfilelines[eachindex].split()#we can select the model that has the closest value for the y co-ordinate when compared to the value stored
								a=eachindex				#in yco
								tempcoords=open(path+'tempcoords.txt','a')
								if records[0]=='MODEL' and float(records[2][:7])==smallest:
									while libfilelines[a].split()[0]!='ENDMODEL':
										if libfilelines[a].split()[0]=='MODEL':
											tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
										a+=1
								tempcoords.close()

								smallesty=999999999
								smallestyindex=-1
								tempcoords=open(path+'tempcoords.txt','r')
								templines=tempcoords.readlines()
								for eachrecord in templines:				#Finds the model number of the model with the closest value for y co-ordinate
									splitrec=eachrecord.split()
									if abs(float(splitrec[1])-yco)<abs(smallesty):
										smallesty=float(splitrec[1])-yco
										smallestyindex=splitrec[0]
								tempcoords.close()

							for eachindex in range(len(libfilelines)):		#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
								records=libfilelines[eachindex].split()		#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
								if records[0]=='MODEL':				#closest to the value in yco
									if float(records[2][:7])==smallest:
										if records[1]==smallestyindex:
											a=eachindex
											while True:
												libfilerecords=libfilelines[a].split()
												if libfilerecords[0]=='MODEL':
													pass
												elif libfilerecords[0]=='ENDMODEL':
													break
												else:
													if libfilerecords[1]=='N':
														a+=1
														while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file
															a+=1
													libfilerecords=libfilelines[a].split()
													to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
													if libfilerecords[1]=='CA':
														xco+=float(libfilerecords[4])
														yco+=float(libfilerecords[5])
														zco+=float(libfilerecords[6])
													newpdb.write('\t'.join(to_write))
													tracker+=1
												a+=1
											break
							libfile.close()
							aminonum+=1
							print '_',
							j+=1
							clusters=[]
						newpdb.close()
						writePDBFormat(path+'temp.pdb',path+'WXYZ.pdb')
						exit()
			aminonum+=1
			print query[j],
			j+=1
			k+=1
			l+=1
			m+=1
		else:
			libfile=open(lib2path+query[j-1]+'_'+query[j]+'.txt','r')
			libfilelines=libfile.readlines()
			clusterline=libfilelines[0]
			smallest=99999999
			for i in range(len(clusterline.split())):
				if i!=0:
					clusters.append(float(clusterline.split()[i]))
			for each in clusters:
				if abs(each-xco)<abs(smallest-xco):
					smallest=float(each)			#identifies the cluster closest to the xco of the last CA

			tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file
			tempcoords.close()

			for eachindex in range(len(libfilelines)):			#From the library file, y co-ordinates of all the models in the closest cluster are copied to a temp file so that we can
				records=libfilelines[eachindex].split()			#select the model that has the closest value for the y co-ordinate when compared to the value stored in yco
				a=eachindex
				tempcoords=open(path+'tempcoords.txt','a')
				if records[0]=='MODEL' and float(records[2][:7])==smallest:
					while libfilelines[a].split()[0]!='ENDMODEL':
						if libfilelines[a].split()[0]=='MODEL':
							tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
						a+=1
				tempcoords.close()

				smallesty=999999999
				smallestyindex=-1
				tempcoords=open(path+'tempcoords.txt','r')
				templines=tempcoords.readlines()
				for eachrecord in templines:				#Finds the model number of the model with the closest value for y co-ordinate
					splitrec=eachrecord.split()
					if abs(float(splitrec[1])-yco)<abs(smallesty):
						smallesty=float(splitrec[1])-yco
						smallestyindex=splitrec[0]
				tempcoords.close()

			for eachindex in range(len(libfilelines)):		#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
				records=libfilelines[eachindex].split()		#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
				if records[0]=='MODEL':				#closest to the value in yco
					if float(records[2][:7])==smallest:
						if records[1]==smallestyindex:
							a=eachindex
							while True:
								libfilerecords=libfilelines[a].split()
								if libfilerecords[0]=='MODEL':
									pass
								elif libfilerecords[0]=='ENDMODEL':
									break
								else:
									if libfilerecords[1]=='N':
										a+=1
										while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file as we have already written it
											a+=1
									libfilerecords=libfilelines[a].split()
									to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
									if libfilerecords[1]=='CA':
										xco+=float(libfilerecords[4])
										yco+=float(libfilerecords[5])
										zco+=float(libfilerecords[6])
									newpdb.write('\t'.join(to_write))
									tracker+=1
								a+=1
							break
			aminonum+=1
			libfile.close()
			clusters=[]
			print '_',
			if blast_record.alignments[0].hsps[0].sbjct[k]=='-':	#skips through the hyphens in the subject sequence of the BLAST result
				j+=1
				k+=1
				l+=1
			else:
				j+=1
				k+=1
				l+=1
				m+=1

	while j<input_fasta_len:
		libfile=open(lib2path+query[j-1]+'_'+query[j]+'.txt','r')
		libfilelines=libfile.readlines()
		clusterline=libfilelines[0]
		smallest=99999999
		for i in range(len(clusterline.split())):
			if i!=0:
				clusters.append(float(clusterline.split()[i]))
		for each in clusters:
			if abs(each-xco)<abs(smallest-xco):
				smallest=float(each)			#identifies the cluster closest to the xco of the last CA

		tempcoords=open(path+'tempcoords.txt','w+')		#creating a temp file
		tempcoords.close()

		for eachindex in range(len(libfilelines)):			#From the library file, y co-ordinates of all the models in the closest cluster are copied to a temp file so that we can
			records=libfilelines[eachindex].split()			#select the model that has the closest value for the y co-ordinate when compared to the value stored in yco
			a=eachindex
			tempcoords=open(path+'tempcoords.txt','a')
			if records[0]=='MODEL' and float(records[2][:7])==smallest:
				while libfilelines[a].split()[0]!='ENDMODEL':
					if libfilelines[a].split()[0]=='MODEL':
						tempcoords.write(libfilelines[a].split()[1]+'\t'+libfilelines[a+2].split()[5]+'\n')
					a+=1
			tempcoords.close()

			smallesty=999999999
			smallestyindex=-1
			tempcoords=open(path+'tempcoords.txt','r')
			templines=tempcoords.readlines()
			for eachrecord in templines:				#Finds the model number of the model with the closest value for y co-ordinate
				splitrec=eachrecord.split()
				if abs(float(splitrec[1])-yco)<abs(smallesty):
					smallesty=float(splitrec[1])-yco
					smallestyindex=splitrec[0]
			tempcoords.close()

		for eachindex in range(len(libfilelines)):		#From the library file, a model is selected to fill the gap in the query match if that model is in the cluster that is closest
			records=libfilelines[eachindex].split()		#to the value in xco and if it's model number is equal to the model number of the model which has a CA with y co-ordinate value
			if records[0]=='MODEL':				#closest to the value in yco
				if float(records[2][:7])==smallest:
					if records[1]==smallestyindex:
						a=eachindex
						while True:
							libfilerecords=libfilelines[a].split()
							if libfilerecords[0]=='MODEL':
								pass
							elif libfilerecords[0]=='ENDMODEL':
								break
							else:
								if libfilerecords[1]=='N':
									a+=1
									while libfilelines[a].split()[1]!='N':		#skipping the first amino acid in the library file as we have already written it
										a+=1
								libfilerecords=libfilelines[a].split()
								to_write=[libfilerecords[0],str(tracker+1),libfilerecords[1],libfilerecords[2],libfilerecords[3],str(aminonum),str(float(libfilerecords[4])+xco),str(float(libfilerecords[5])+yco),str(float(libfilerecords[6])+zco),libfilerecords[7],libfilerecords[8],libfilerecords[9],'\n']
								if libfilerecords[1]=='CA':
									xco+=float(libfilerecords[4])
									yco+=float(libfilerecords[5])
									zco+=float(libfilerecords[6])
								newpdb.write('\t'.join(to_write))
								tracker+=1
							a+=1
						break
		aminonum+=1
		libfile.close()
		clusters=[]
		print '_',
		j+=1
usepdb.close()
newpdb.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Writing in pdb format

writePDBFormat(path+'temp.pdb',path+'WXYZ.pdb')

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
