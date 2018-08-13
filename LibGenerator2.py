import os
import math
import glob
import amino
import numpy as np
from sklearn.cluster import KMeans

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#initializing variables

path='/home/intern/Desktop/UpdatedPDB/'
libpath='/home/intern/Desktop/Library2/'
L=[]
A=amino.aminotuple(2)
modelno=1
tempinput=[]
clusterinput=[]
xco=yco=zco=0.0
smallest=9999999999

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for eachamino in A:							#iterates for each of the 400 combinations
	print eachamino,'\n'
	newfile=open(libpath+'_'.join(eachamino)+'.txt','w+')
	modelno=1
	for filename in glob.iglob(os.path.join(path,'*.pdb')):		#iterates through each pdb file in the folder for every amino acid
		L=[]				
		handle=open(filename,'r')
		lines=handle.readlines()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		for i in range(len(lines)-1):
			handle=lines[i]
			currentline=handle.split()			#storing indices of all the alpha carbons in a list
			types=currentline[2]				
			title=currentline[0]				
			if types=='CA':
				L.append(i)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		for i in range(len(L)-1):			#iterates for all the alpha carbons
			t1=L[i]
			t2=L[i+1]
			current,next=lines[t1],lines[t2]
			if current.split()[3]==eachamino[0] and next.split()[3]==eachamino[1]:	#checks equality of pair of amino acids
				itervar=t1-1
				iterhandle=lines[itervar]
				iterator=iterhandle.split()
				chainid=iterator[4]

				xco=float(current.split()[6])			#\
				yco=float(current.split()[7])			# - storing the x,y,z coordinates of the alpha carbon
				zco=float(current.split()[8])			#/

				model=['MODEL',str(modelno),'\n']
				newfile.writelines(' '.join(model))
				flag=0

				while flag!=2:					#iterates till the 3rd N is encountered, where the first 2 Ns belong to each of the amino acids in the pair we have taken
										#this is done because for the alpha carbon, the coordinates should be written as it is and relative distance need not be
					if itervar==t1:				#calculated, as clustering is based on the x,y,z coordinates of the alpha carbon
						relativex=xco
						relativey=yco
						relativez=zco

					else:
						relativex=float(iterator[6])-xco	#\
						relativey=float(iterator[7])-yco	# -calculating the relative distance from CA
						relativez=float(iterator[8])-zco	#/

					try:
						to_write=[iterator[0],iterator[2],iterator[3],iterator[4],str(relativex),str(relativey),str(relativez),iterator[9],iterator[10],iterator[11],iterator[12],'\n']

					except:
						to_write=[iterator[0],iterator[2],iterator[3],iterator[4],str(relativex),str(relativey),str(relativez),iterator[9],iterator[10],iterator[11],'\n']
					newfile.writelines(' '.join(to_write))
					itervar+=1

					try:
						iterhandle=lines[itervar]
						iterator=iterhandle.split()		#moving to the next line

					except IndexError:
						break

					if iterator[0]=='ATOM' and iterator[2]=='N':		#flag increased by 1 when next N is encountered
						flag+=1

				newfile.write('ENDMODEL\n')
				modelno+=1
	newfile.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#KMeans clustering

for filename in glob.iglob(os.path.join(libpath,'*.txt')):			#takes each file in the library
	print filename
	modelno=1
	tempinput=[]
	clusterinput=[]
	handle=open(filename,'r')
	lines=handle.readlines()
	for i in range(len(lines)):
		modelno=1
		if lines[i].split()[0]=='MODEL':
			tempinput.append(float(lines[i+2].split()[4]))		#extracting the x coordinate of all the alpha carbons
			clusterinput.append(tempinput)
			tempinput=[]
	handle.close()
	kmeans=KMeans(n_clusters=7)

	try:
		kmeans.fit(np.array(clusterinput))				#clustering

	except:
		continue

	centroids=kmeans.cluster_centers_					#returns the centroids

	temphandle=open(libpath+'temp.txt','w+')
	to_write=['CLUSTER',str(centroids[0][0])[:7],str(centroids[1][0])[:7],str(centroids[2][0])[:7],str(centroids[3][0])[:7],str(centroids[4][0])[:7],str(centroids[5][0])[:7],str(centroids[6][0])[:7],'\n']	#writing the centroids at the start of the file
	temphandle.writelines('\t'.join(to_write))

	for eachcluster in centroids:
		to_write=['CLUSTER',str(eachcluster[0]),'\n']
		temphandle.writelines('\t'.join(to_write))

		for i in range(len(lines)):
			if lines[i].split()[0]=='MODEL':
				smallest=9999999999999999
				xco=float(lines[i+2].split()[4])
				for j in range(7):
					if abs(xco-centroids[j][0])<=abs(xco-smallest):		#identifying which cluster the model belongs to
						smallest=centroids[j][0]

				if smallest==eachcluster[0]:
					to_write=['MODEL',str(modelno),str(eachcluster[0]),'\n']
					temphandle.write(' '.join(to_write))
					i+=1
					
					while lines[i].split()[0]!='ENDMODEL':
						temphandle.write(lines[i])
						i+=1

					to_write='ENDMODEL\n'
					temphandle.write(to_write)
					modelno+=1
		temphandle.write('ENDCLUSTER\n')

	temphandle.close()
	os.rename(libpath+'temp.txt',filename)				#replacing the original file with a file that is clustered

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
