import networkx as nx
import numpy as np
from sklearn.cluster import SpectralClustering
import networkx.algorithms as nxAlg
from numpy import linalg as LA
import math
from scipy.sparse import csgraph
from numpy import genfromtxt
import pandas as pd
import matplotlib.pyplot as plt
import os

def CheckConnectiveness(Ajacency):
	FlagConnect=False
	Laplacinan=csgraph.laplacian(Ajacency, normed=False)
	w_eig,v_eig=LA.eig(Laplacinan)
	eigenLaplacian=np.absolute(w_eig)
	eigenLaplacian.sort()
	lambda1=eigenLaplacian[1]	
	NotconnectedGRaph=lambda1<1e-15
	if NotconnectedGRaph:
		FlagConnect=False
		#print("Graph is Not Connected" )
	else:
		FlagConnect=True
		#print("GRaph IS connected :) ")
	return FlagConnect,lambda1
	
def RamanujanCheck(Adjacency,degree):
	FlagRam=False
	w,v=LA.eig(Adjacency)
	eigenAdj=np.absolute(w)
	eigenAdj.sort()
	sizeVec=len(eigenAdj)

	mu0=eigenAdj[N-1]
	mu1=eigenAdj[N-2]

	RamanujanVal=2*(math.sqrt(degree-1))
	RamanujanVal=2*(math.sqrt(degree-1))
	Ramanujan=mu1 <= RamanujanVal
	#print(mu1, "<=", RamanujanVal)
	if (Ramanujan):
		FlagRam=True
	else:
		FlagRam=False
	return FlagRam	
	
def Ncut_conduct_Graph(Graph,Ajacency,name):
	clustering1=SpectralClustering(2,affinity='precomputed', n_init=1000)
	clustering1.fit(Ajacency)
	clusterLabels1=clustering1.labels_
	cluster1=[index for index in range(len(clusterLabels1)) if clusterLabels1[index] == 0]
	cluster2=[index for index in range(len(clusterLabels1)) if clusterLabels1[index] == 1]
	subset1=len(cluster1)
	subset2=len(cluster2)
	#print("cluster1_",name, ": ",cluster1,"cluster2_",name,": ",cluster2)
	Vol_1=nxAlg.volume(Graph,cluster1)
	Vol_2=nxAlg.volume(Graph,cluster2)
	cutSize=nxAlg.cut_size(Graph,cluster1,cluster2)
	conductance=nxAlg.conductance(Graph,cluster1,cluster2)
	Ncut=nxAlg.normalized_cut_size(Graph,cluster1,cluster2,weight=None)
	return Vol_1,Vol_2,cutSize,Ncut,conductance
	
def GraphCreation(type,Nodes,Edges,degree):
	ZeroVal=True
	NotConnectFlag=True
	RamFlag=False
	iter=0
	while ZeroVal or  NotConnectFlag or (not RamFlag):
		iter=iter+1
		print("--------iter",iter,"-------")
		#print("Zero in Degree list: ",ZeroVal, ", Connectivity",NotConnectFlag)
		if (type=="Regular"):
			Graph=nx.random_regular_graph(degree,Nodes)
		elif (type=="Erdos"):
			RamFlag=True
			Graph=nx.gnm_random_graph(N,edges)
		else:
			print("!!!!!!!Choose correct Graph Type!!!!!!!")
		Ajacency=(nx.to_numpy_matrix(Graph)).astype(int)
		Ajacency=np.matrix(Ajacency)
		DegreeList=sum(np.array(sum(Ajacency)))
		ZeroVal=0 in DegreeList
		ConnectFlag,lambda_1 =CheckConnectiveness(Ajacency)
		NotConnectFlag=ConnectFlag==False
		if (type=="Regular"):
			RamFlag=RamanujanCheck(Ajacency,sum(DegreeList)/len(Ajacency))
	return Ajacency, DegreeList,Graph
	


def CheckSymmetric(Adjacency):
	Symetric=np.allclose(Adjacency,Adjacency.T,atol=1e-15)
	return Symetric
	
def ReadGraphGiveConnectivity(iGraphType, iNodes, iFolder,iDegree):
	File_path=	Folder_path+ str(iNodes)+"/"+iGraphType +str(iFolder)+ "/"+str(iDegree)+iGraphType+".graph"

	print("------ ", G,' ,folder ', str(iFolder), " ,degree ",str(iDegree), "-------- ")		
	A=genfromtxt(File_path)
	Ajacency=np.array(A)
	Connect_Flag,lambda1=CheckConnectiveness(Ajacency)
	SymetricMat=CheckSymmetric(Ajacency)
	if (Connect_Flag and SymetricMat):
		print("connectivity is: ",lambda1 )
		dataframe = pd.DataFrame(data=Ajacency)
		degrees=sum(np.array(sum(Ajacency)))
		degree=degrees.mean().astype(int)
	else:
		print('******************** Error: Graph does not have our criteria**************************')
	return lambda1

def SaveFile(iDirectory,iName):
    print(iDirectory)
    if not os.path.exists(iDirectory):
    	os.makedirs(iDirectory)
    MyFormat=['eps','pdf']
    for form in MyFormat:
    	savingPath=iDirectory+iName+'.'+ form
    	manager=plt.get_current_fig_manager()
    	manager.window.showMaximized()
    	plt.savefig(savingPath, bbox_inches = "tight", dpi=500)
################# Inputs######################################
NoNodes=120
if NoNodes==240:
	Degree=[90,100,110,120,130,140,150,160,170,180,190,200]
	myXlabels=['90','100','110','120','130','140','150','160','170','180','190','200','239']
elif NoNodes==40:
	Degree=[7,11,15,19,21,23,27,31,39]
	myXlabels=['7','11','15','19','21','23','27','31','39']
elif NoNodes==120:
	Degree=[10,20,30,40,50,60,70,80,90,100]
	myXlabels=['10','20','30','40','50','60','70','80','90','100','119']
else:
	print('Wrong graphs')


Folder_path="./../graphs/"
GraphType=['Regular', 'Mean']
fontSize=16
i=0
plt.figure()
#plt.title('Connectivity')
marker='*'
color='r'
myplot=[]
labels=[]
complete=NoNodes-1
folder=100 #number of batch of graphs that i want to plot, here i want to plot 100 regular and mean graphs
for G in GraphType:
	PlotName='Connectivity for '+G+' graphs'
	labels.append(G)
	for F in range(1,folder+1):
		connectivity=[]
		for D in Degree:
			lambda_1=ReadGraphGiveConnectivity(G,NoNodes, F, D )
			connectivity.append(lambda_1)
		if i==0 or i ==folder: # size of i should be in the size of folders for regular or mean graphs and not for degrees
			if G=='Regular':
				input_marker='o'
			else:
				input_marker='*'
			print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" , i , " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
			plt.scatter(Degree, connectivity, facecolors='none', label=G ,edgecolors=color, marker=input_marker)
		else:
			if G=='Regular':
				marker='o'
			else:
				marker='*'
			plt.scatter(Degree, connectivity, facecolors='none', edgecolors=color, marker=input_marker)
		i=i+1

			
	marker='o'
	color='b'
lamda_complete=ReadGraphGiveConnectivity('Regular',NoNodes, 1, complete )
plt.scatter(complete, lamda_complete, facecolors='none', label="Complete", edgecolors='k', marker='^')		
#plt.title('Connectivity')
plt.ylabel('Connectivity',{"size":fontSize})
plt.xlabel('Degree',{"size":fontSize-5})
Degree.append(complete)
print('degree ' ,Degree)
print('label ',myXlabels)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(Degree,myXlabels)
plt.tick_params(labelsize=fontSize-5)
plt.legend(prop={"size":fontSize})
myDirectory='./../pictures/'
name='ConnectivityOf100Reg_Mean_CompleteFor'+str(NoNodes)+'Nodes'
SaveFile(myDirectory,name)
plt.show()