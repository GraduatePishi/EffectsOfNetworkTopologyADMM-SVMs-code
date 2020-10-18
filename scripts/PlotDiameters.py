# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:14:09 2020

@author: SHTA
"""
import matplotlib.pyplot as plt
import numpy as np 
from pathlib import Path

def ReadMyFile(iFileName):
    myFile=Path(iFileName)
    if not myFile.exists(): 
        print(iFileName , 'does NOT exist!!!')
    token = open(iFileName,'r')
    linestoken=token.readlines()[1:]

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

path='./../graphs/240/diameter/5000Diameters/Dimater_for_'
Degree=[3,6,12,24,48,96]

r=0
pos=[]
while r < len(Degree):
    pos.append(r)
    r=r+1
fontSize  =16 
TotalDiameter=[]
TotalLambda=[]
for count, D in enumerate(Degree):
    Diameter=[]
    Lambdas=[]
    pathTotal=path+str(D)+'-degree'
    Diameter,Lambdas=ReadMyFile(pathTotal)
    print('-----------',str(D),' Degree---------')
    print('Smallest Diameter:', min(Diameter), ', Largest Diameter:', max(Diameter))
    print('Smallest Lambdas:', min(Lambdas), ', Largest Lambdas:', max(Lambdas))
    print()
    TotalDiameter.append(Diameter)
    TotalLambda.append(Lambdas)
    
#plt.title('The Diameters of mean degree graph',{"size":fontSize+2})   
box=plt.boxplot(TotalDiameter, positions=pos) 
plt.xlabel('Degree',{"size":fontSize+2})
plt.ylabel('Diameter',{"size":fontSize+2})
plt.xticks(pos, Degree)
plt.tick_params(labelsize=fontSize)
manager=plt.get_current_fig_manager()
manager.window.showMaximized()
savingPath='./../pictures/Diameter.eps'
plt.savefig(savingPath, bbox_inches = "tight", dpi=500)
plt.show()

#plt.title('The Connectivity of mean degree graphs',{"size":fontSize+2})
box=plt.boxplot(TotalLambda, positions=pos) 
plt.xlabel('Degree',{"size":fontSize+2})
plt.ylabel('Connectivity',{"size":fontSize+2})
plt.xticks(pos, Degree)
plt.tick_params(labelsize=fontSize)
manager=plt.get_current_fig_manager()
manager.window.showMaximized()
savingPath='./../pictures/Connectivity.eps'
plt.savefig(savingPath, bbox_inches = "tight", dpi=500)
plt.show()