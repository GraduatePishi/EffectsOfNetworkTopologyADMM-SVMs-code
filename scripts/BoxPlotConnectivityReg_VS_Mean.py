# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:54:11 2020

@author: SHTA
"""


import numpy as np
from numpy import linalg as LA
from scipy.sparse import csgraph
from numpy import genfromtxt
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path

def ReadMyFile(iFileName):
    myFile=Path(iFileName)
    if not myFile.exists(): 
        print(iFileName , 'does NOT exist!!!')
    token = open(iFileName,'r')
    linestoken=token.readlines()

    FirstCol=[]
    SecondCol=[]
    NrOdCols=2
    for x in linestoken:
        value=x.split()
        elements=len(x.split())
        for i in range(0,elements):
            if i%NrOdCols ==0:
                FirstCol.append(int(value[i]))
            else:
                SecondCol.append(float(value[i]))
    token.close()
    
    return FirstCol,SecondCol

def CheckConnectiveness(Ajacency):
    FlagConnect=False
    Laplacinan=csgraph.laplacian(Ajacency, normed=True)
    w_eig,v_eig=LA.eig(Laplacinan)
    eigenLaplacian=np.absolute(w_eig)
    eigenLaplacian.sort()
    lambda1=eigenLaplacian[1]    
    NotconnectedGRaph=lambda1<1e-15
    if NotconnectedGRaph:
        FlagConnect=False
        print("Graph is Not Connected" )
    else:
        FlagConnect=True
    return FlagConnect,lambda1
    
def CheckSymmetric(Adjacency):
    Symetric=np.allclose(Adjacency,Adjacency.T,atol=1e-15)
    return Symetric
    
def ReadGraphGiveConnectivity(iGraphType, iNodes, iFolder,iDegree):
    File_path=    Folder_path+iGraphType +str(iFolder)+ "/"+str(iDegree)+iGraphType+".graph"
    A=genfromtxt(File_path)
    Ajacency=np.array(A)
    Connect_Flag,lambda1=CheckConnectiveness(Ajacency)
    SymetricMat=CheckSymmetric(Ajacency)
    if (Connect_Flag and SymetricMat):
        dataframe = pd.DataFrame(data=Ajacency)
        degrees=sum(np.array(sum(Ajacency)))
        degree=degrees.mean().astype(int)
    else:
        print('******************** Error: Graph does not have our criteria**************************')
    return lambda1
################# Inputs######################################
fontSize=16
NoNodes=240

AllDegrees=[3,6,12,24,48,96]


Folder_path="./../graphs/"+str(NoNodes) +"/Connectivity_Comparison_Reg_vs_Mean/"

GraphType=['Regular', 'Mean']

i=0
plt.figure()
marker='*'
color='r'
myplot=[]
labels=[]
myXlabels=[]

complete=NoNodes-1

folder=100#number of batch of graphs that i want to plot, here i want to plot 100 regular and mean graphs
for Degree in AllDegrees:
    degrees_reg=[ i for i in range(1, 1+1)]
    ones=np.ones(len(degrees_reg))*0.3
    degrees_mean=degrees_reg+ones
    for G in GraphType:
        Mycolor=[]
        Mylabel=[]
        PlotName='Connectivity for '+G+' graphs'
        labels.append(G)
        if G=='Regular':
            MaxConnectivityForDegree=np.zeros(1)
        else:
            MaxConnectivityForDegree=np.ones(1)
        FolderNameOfMaxConnect=np.zeros(1)
        SaveAllconnectivity=[]

        connectivity=[]
        myXlabels=[]
        myXlabels.append(str(Degree))
        i=0
        if G=='Regular':
            for F in range(1,folder+1):
                lambda_1=ReadGraphGiveConnectivity(G,NoNodes, F, Degree )
                connectivity=connectivity+[lambda_1]
        else:
            print('connectivity for degree ', Degree,' and ', G, 'in folder', F)#, ' is: ',connectivity)
            path='./../graphs/240/diameter/Dimater_for_'
            pathTotal=path+str(Degree)+'-degree'
            Diameter,Lambdas=ReadMyFile(pathTotal)
            connectivity=Lambdas

        SaveAllconnectivity.append(connectivity)
        Mylabel.append(str(Degree))
        
        if G=='Regular':
            colorCode='-r*'
            x=Degree    
            pos=degrees_reg
            Mycolor.append('white')
        elif G=='Mean':
            colorCode='-b*'
            x=Degree+0.3
            pos=degrees_mean
            Mycolor.append('black')
        else:
            print("WRONG color code")
                
    
        if G=='Regular':
            RegBox1=plt.boxplot(SaveAllconnectivity, positions=pos, patch_artist=True, labels= Mylabel)
            for patch, color in zip(RegBox1['boxes'], Mycolor):
                patch.set_facecolor(color)
        elif G=='Mean':
            RegBox2=plt.boxplot(SaveAllconnectivity, positions=pos, patch_artist=True, labels= Mylabel)
            for patch, color in zip(RegBox2['boxes'], Mycolor):
                patch.set_facecolor(color)
    
    
        marker='o'
        color='b'
    CompeleteFlag=False
    if CompeleteFlag==True:
        lamda_complete=ReadGraphGiveConnectivity('Regular',NoNodes, 1, complete )
        plt.scatter(complete, lamda_complete, facecolors='none', label="Complete", edgecolors='k', marker='^')    
        Degree.append(complete)
    
    plt.ylabel('Connectivity',{"size":fontSize})
    plt.xlabel('Degree',{"size":fontSize-2})
    
    print('degree ' ,Degree)
    
    plt.legend([RegBox1["boxes"][0], RegBox2["boxes"][0]], GraphType, prop={"size":fontSize})#,loc='upper left'
    ax=plt.gca() #for creating the defined x-axes
    plt.tick_params(labelsize=fontSize)
    
    MyFormat=['eps','pdf']
    myDirectory='./../pictures/'
    if not os.path.exists(myDirectory):
        os.makedirs(myDirectory)
    for form in MyFormat:
            savingPath=myDirectory+'BoxConnectivity_'+str(NoNodes)+'N_'+str(Degree)+'.'+ form
            manager=plt.get_current_fig_manager()
            manager.window.showMaximized()
            plt.savefig(savingPath, bbox_inches = "tight", dpi=500)
    
    plt.show()
