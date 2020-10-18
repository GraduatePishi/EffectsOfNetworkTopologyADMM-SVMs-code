
"""
@author: Shirin Tavara
"""

import matplotlib.pyplot as plt
import os
from pathlib import Path


def ReadMyFile(iFileName,minIter,myIterNumber):
    myFile=Path(iFileName)
    if not myFile.exists(): 
        print(iFileName , 'does NOT exist!!!')
    token = open(iFileName,'r')
    linestoken=token.readlines()[1:]
    NrOdCols=12
    iter=[]
    Wp=[]
    Bp=[]
    Wd=[]    
    Bd=[]
    maxVall_Vj=[]
    maxVi_Vj=[]
    AVGVall_Vj=[]
    AVGVi_Vj=[]
    AVGTime=[]
    MAXsolverIter=[]
    SolvTime=[]
    
#Order of variables: Iter Wp Bp Wd Bd max(Vall-Vj) max(Vj-Vi) Avg(Vall-Vj) Avg(Vj-Vi) Avg(time) MaxSolvIter SolverTime
    for x in linestoken:
        value=x.split()
        elements=len(x.split())
        for i in range(0,elements):
            if i%NrOdCols ==0:
                iter.append(int(value[i]))
            elif i%NrOdCols==1:
                Wp.append(float(value[i]))
            elif i%NrOdCols==2:
                Bp.append(float(value[i]))
            elif i%NrOdCols==3:
                Wd.append(float(value[i]))
            elif i%NrOdCols==4:
                Bd.append(float(value[i]))
            elif i%NrOdCols==5:
                maxVall_Vj.append(float(value[i]))
            elif i%NrOdCols==6:
                maxVi_Vj.append(float(value[i]))
            elif i%NrOdCols==7:
                AVGVall_Vj.append(float(value[i]))
            elif i%NrOdCols==8:
                AVGVi_Vj.append(float(value[i]))
            elif i%NrOdCols==9:
                AVGTime.append(float(value[i]))
            elif i%NrOdCols==10:
                MAXsolverIter.append(float(value[i]))
            elif i%NrOdCols==11:
                SolvTime.append(float(value[i]))

    token.close()
    print('MaxIteration',myIterNumber)
    return iter[minIter:myIterNumber],Wp[minIter:myIterNumber],Bp[minIter:myIterNumber],Wd[minIter:myIterNumber], Bd[minIter:myIterNumber], maxVall_Vj[minIter:myIterNumber], maxVi_Vj[minIter:myIterNumber], AVGVall_Vj[minIter:myIterNumber],AVGVi_Vj[minIter:myIterNumber], AVGTime[minIter:myIterNumber], MAXsolverIter[minIter:myIterNumber] ,SolvTime[minIter:myIterNumber]


def SaveFile(iDirectory,iName):
    print(iDirectory)
    if not os.path.exists(iDirectory):
        os.makedirs(iDirectory)
    MyFormat=['eps','pdf']
    for form in MyFormat:
        savingPath=iDirectory+iName+'.'+ form
        manager=plt.get_current_fig_manager()
        manager.window.showMaximized()
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        plt.savefig(savingPath, bbox_inches = "tight", dpi=450)
    
NoNodes=240
fontSize=16
MinIter=2
MaxIter=24

Degree=[24,48,96]#[3,4,6,12]#
MyPlot=[3,4]
DatasetName='higgs1m' 
#################### READ ME ######################
if DatasetName=='higgs1m':
    if Degree[0]==3:
        MinIter=0
        MaxIter=18
    else:
        MinIter=0
        MaxIter=18       
    
elif DatasetName=='susy700k':
    MinIter=0
    MaxIter=5
    
elif DatasetName=='covtype400k':
    MinIter=0
    MaxIter=12
myColor=['r', 'b', 'k', 'g', 'y', 'm','c','brown']
colorIndex=0
myXlabels=[]
mak=['*', 'o', 'v','s' ,'D','8']

Total_max_Vall_Vj=[] #1
Total_max_Vj_Vi=[]#2
Total_Avg_Vall_Vj=[]#3
Total_Avg_Vj_Vi=[]#4
Total_Avg_time=[]#5
Total_MaxSolvIter=[]#6
Total_SolvTime=[]#7
fontSize=16
MyDash=['solid','dashed','dashdot','solid', 'dotted','dashed','dashdot']
Mylegend=[]
myPlot=[]


filePath="./../Results/"+ DatasetName+"/"

constant=1
Length=0
GraphType=['Regular', 'Mean']
if Degree[0]< 20:
    OtherGraphs=['24_10dim_2DTorus','Ring','Line']
else:
    OtherGraphs=[]
    


for plotCounter, plott in enumerate(MyPlot):
    MyTitle1="_"
    myCounter=0
    for counter, D in enumerate(Degree):
        myCounter=counter
        Total_max_Vall_Vj=[] #1
        Total_max_Vj_Vi=[]#2
        Total_Avg_Vall_Vj=[]#3
        Total_Avg_Vj_Vi=[]#4
        Total_Avg_time=[]#5
        Total_MaxSolvIter=[]#6
        Total_SolvTime=[]#7
        Total_Iter=[]
        for count, G in enumerate(GraphType):
            print("Dataset is ", DatasetName, " ----- Node ", NoNodes, " -----", D, G )
            fileName=filePath+str(D)+G+"_"+str(NoNodes)+"N/All Max and AVG Erros for " +str(D)+ G+".error2"

            iteration=[]
            W_prim=[]
            B_prim=[]
            W_dual=[]
            B_dual=[]
            max_Vall_Vj=[] #1
            max_Vj_Vi=[]#2
            Avg_Vall_Vj=[]#3
            Avg_Vj_Vi=[]#4
            Avg_time=[]#5
            MaxSolvIter=[]#6
            SolvTime=[]#7
    
            iteration, W_prim,B_prim, W_dual, B_dual, max_Vall_Vj,max_Vj_Vi,Avg_Vall_Vj, Avg_Vj_Vi, Avg_time,MaxSolvIter, SolvTime =ReadMyFile(fileName,MinIter,MaxIter)
            
            Total_max_Vall_Vj.append(max_Vall_Vj)
            Total_max_Vj_Vi.append(max_Vj_Vi)
            Total_Avg_Vall_Vj.append(Avg_Vall_Vj)
            Total_Avg_Vj_Vi.append(Avg_Vj_Vi)
            Total_Avg_time.append(Avg_time)
            Total_MaxSolvIter.append(MaxSolvIter)
            Total_SolvTime.append(SolvTime)
            Total_Iter.append(iteration)
    
           
            print('----------- DONE ------------')
            print()

            Mylegend.append(str(D)+G)
            if G=='Regular':
                MyTitle1=MyTitle1+ str(D)+'_' 
            

        print(MyTitle1)

        myDirectory='./../pictures/'+DatasetName+'/'
    
            
        plt.rc('font', size=fontSize-4)
        #print('iteration',iteration)
        NrGraphsTypes=len(GraphType)
        Mylabel=GraphType
        if OtherGraphs !=[] and D==Degree[-1]:
            for Othercount, O in enumerate(OtherGraphs):
                NrGraphsTypes=NrGraphsTypes+1
                print('I AM HEREE!!!!!!!!!!!!!')
                print("Dataset is ", DatasetName, " ----- Node ", NoNodes, " -----", OtherGraphs[Othercount] )
                fileName=filePath+OtherGraphs[Othercount]+"_"+str(NoNodes)+"N/All Max and AVG Erros for " + OtherGraphs[Othercount]+".error2"
                iteration=[]
                W_prim=[]
                B_prim=[]
                W_dual=[]
                B_dual=[]
                max_Vall_Vj=[] #1
                max_Vj_Vi=[]#2
                Avg_Vall_Vj=[]#3
                Avg_Vj_Vi=[]#4
                Avg_time=[]#5
                MaxSolvIter=[]#6
                SolvTime=[]#7
        
                iteration, W_prim,B_prim, W_dual, B_dual, max_Vall_Vj,max_Vj_Vi,Avg_Vall_Vj, Avg_Vj_Vi, Avg_time,MaxSolvIter, SolvTime =ReadMyFile(fileName,MinIter,MaxIter)
                
                Total_max_Vall_Vj.append(max_Vall_Vj)
                Total_max_Vj_Vi.append(max_Vj_Vi)
                Total_Avg_Vall_Vj.append(Avg_Vall_Vj)
                Total_Avg_Vj_Vi.append(Avg_Vj_Vi)
                Total_Avg_time.append(Avg_time)
                Total_MaxSolvIter.append(MaxSolvIter)
                Total_SolvTime.append(SolvTime)
                Total_Iter.append(iteration)
                Mylabel=Mylabel+OtherGraphs
                myCounter=counter+1
                if O=='24_10dim_2DTorus':
                    Mylegend.append('2DTorus')
                else:
                    Mylegend.append(OtherGraphs[Othercount])
            
        if plott==1:
            name='Max_Solv_Iter_lowDegree_'+DatasetName+MyTitle1
            #plt.title('Max iteration of the inner solver',{"size":fontSize})
            for i in range(0,NrGraphsTypes):
                print('Iteration for ', str(i),'th graph, degree ', D, ': ', Total_Iter[i])
                print('MaxSolvIter for ', str(i),'th graph, degree ', D, ': ', Total_MaxSolvIter[i])
                if i>=2:
                    plt.plot(Total_Iter[i], Total_MaxSolvIter[i],color=myColor[counter+i], linestyle=MyDash[i], marker=mak[i], label=Mylabel[i], mfc='none')#6
                else:
                    plt.plot(Total_Iter[i], Total_MaxSolvIter[i],color=myColor[counter], linestyle=MyDash[i], marker=mak[i], label=Mylabel[i], mfc='none')#6
            plt.xlabel('ADMM Iteration',{"size":fontSize})
            plt.ylabel('Solver Iteration',{"size":fontSize})
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.tick_params(labelsize=fontSize)
    
    
        elif plott==2:
            name='SolvTime_lowDegree_'+DatasetName+MyTitle1
            #plt.title('Total training time  ',{"size":fontSize})
            for i in range(0,NrGraphsTypes):
                print('Total_SolvTime',Total_SolvTime[i])
                if i>=2:
                    plt.plot(Total_Iter[i], Total_SolvTime[i],color=myColor[counter+i], linestyle=MyDash[i], marker=mak[i], label=Mylabel[i], mfc='none')#6
                else:
                    plt.plot(Total_Iter[i], Total_SolvTime[i], color=myColor[counter], linestyle=MyDash[i],  marker=mak[i],label=Mylabel[i], mfc='none')#7
            plt.xlabel('ADMM Iteration',{"size":fontSize})
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 
            plt.tick_params(labelsize=fontSize)
            plt.ylabel('Time in sec',{"size":fontSize})

    
        elif plott==3:
            MinIterPlot=0
            if DatasetName=='covtype400k' or DatasetName=='higgs1m':
                MinIterPlot=2
            name='Time_Error_lowDegree_'+DatasetName+MyTitle1
            #plt.title('Difference between the result of node j with neighbors',{"size":fontSize})
            for i in range(0,NrGraphsTypes):
                if i>=2:
                    plt.plot(Total_SolvTime[i][MinIterPlot:], Total_Avg_Vj_Vi[i][MinIterPlot:],color=myColor[counter+i], linestyle=MyDash[i],  marker=mak[i],label=Mylabel[i], mfc='none')#7
                else:
                    plt.plot(Total_SolvTime[i][MinIterPlot:], Total_Avg_Vj_Vi[i][MinIterPlot:],color=myColor[counter], linestyle=MyDash[i],  marker=mak[i],label=Mylabel[i], mfc='none')#7
            plt.xlabel('Time in sec',{"size":fontSize})
            plt.ylabel('e\N{SUBSCRIPT One}',{"size":fontSize})
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.tick_params(labelsize=fontSize)
            plt.legend(Mylegend, prop={"size":fontSize-2})
            
        elif plott==4:
            MinIterPlot=0
            if DatasetName=='covtype400k'or DatasetName=='higgs1m':
                MinIterPlot=2
            
            name='Time_Error_All_Vj_'+DatasetName+MyTitle1
            #plt.title('Difference between the result of node j with the average',{"size":fontSize})
            for i in range(0,NrGraphsTypes):
                if i>=2:
                    plt.plot(Total_SolvTime[i][MinIterPlot:], Total_Avg_Vall_Vj[i][MinIterPlot:],color=myColor[counter+i], linestyle=MyDash[i],  marker=mak[i],label=Mylabel[i], mfc='none')#7
                else:
                    plt.plot(Total_SolvTime[i][MinIterPlot:], Total_Avg_Vall_Vj[i][MinIterPlot:],color=myColor[counter], linestyle=MyDash[i],  marker=mak[i],label=Mylabel[i], mfc='none')#7
            plt.xlabel('Time in sec',{"size":fontSize})
            plt.ylabel('e\N{SUBSCRIPT TWO}',{"size":fontSize})
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.tick_params(labelsize=fontSize)
            plt.legend(Mylegend, prop={"size":fontSize-2})
    
    
    SaveFile(myDirectory,name)
    plt.show() 
    plt.close()
