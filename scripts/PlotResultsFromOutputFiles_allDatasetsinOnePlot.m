clear all
clc
DatasetName={'higgs1m' 'skin_nonskin' 'susy700k' 'covtype400k'};
%DatasetName={'cod-rna'};
r=1;
colors={'b' 'r' 'g' 'k' 'c' };
legendOffFig=[""];   
fontSize=30;
axesSize=15;
titleSize=40;
for data =1:length(DatasetName)
    if strcmp(DatasetName{data},'higgs1m') || strcmp(DatasetName{data},'skin_nonskin') || strcmp(DatasetName{data},'cod-rna') || strcmp(DatasetName{data},'covtype400k')
        G2=[ 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,210,239];
    elseif strcmp(DatasetName{data},'susy700k')
        G2=[110, 120, 130, 140, 150, 160, 170, 180, 190, 200,210,239];
    elseif strcmp(DatasetName{data},'seismic1') || strcmp(DatasetName{data},'miRNAs')
        G2=[11,15,19,23,27,31,39];
    elseif strcmp(DatasetName{data},'seismic')
        G2=[20,30,40,50,60,70,80,90,100,110,119];
    elseif strcmp(DatasetName{data},'covtype400k11')
        G2=[100,150,160,170,180,190,200,210,239];
    else
        disp("Wrong Degrees!!!");
    end
    
    if  strcmp(DatasetName{data},'skin_nonskin')
        legendOffFig(end+1)=strcat('skin-Nonskin');
    else
        legendOffFig(end+1)=strcat(DatasetName{data});
    end
    sizeG=size(G2,2);
    GraphTopol={'Regular'};
    NumOfFolder=1;
    FixedPath=fullfile('Z:\ShirinsResearchResults\PeerJ\','Results','onlyOutputResults');
    for type=1:length(GraphTopol)
        if type==1
            FolderName={'Regular1'};
            Iter_Regular=[];
        elseif type==2
            FolderName={'Mean1'};
            Iter_mean=[];
        else
            disp('Wrong Graph type!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        end
        
        Separator_=['*****************', num2str(DatasetName{data}) ,'***********************'];
        disp(Separator_);
        iter_tmp_phish=sizeG;
        time_tmp_phish=sizeG;
        
        for k=1:length(FolderName)
            shuffle_svmguid=['*****************',FolderName{k},' ***********************'];
            disp(shuffle_svmguid);
            graphTypeFolder=[GraphTopol{type},num2str(NumOfFolder)];
            filePath=[FixedPath,'\',FolderName{k},'_',DatasetName{data},'\'];
            iterReg=[];
            MasterTimeFromMain=[];
            TotalTrainTimeFromMain=[];
            TotalSolvTimeFromMain=[];
            TotalSolvTimeFromSolver=[];
            NLOPTTime=[];
            CommuniTime=[];
            DistDataTime=[];
            timeMisc=[];
            testing=[];
            scaledCommuniTime=[];
            
            for i=1:length(G2)
                fileOutput1=[FolderName{k},'_', num2str(DatasetName{data}),'*_',num2str(G2(i)),GraphTopol{type},'_*'];
                outputDegree=['... Reading File ', fileOutput1 ];
                disp(outputDegree);
                MyList=dir(fullfile(filePath, fileOutput1));
                MyList([MyList.isdir])=[];
                filePhish=MyList.name;
                cd(filePath);
             
                it=importdata(filePhish);
                iterReg(end+1)=it(1);
                MasterTimeFromMain(end+1)=it(2);
                TotalTrainTimeFromMain(end+1)=it(3);
                TotalSolvTimeFromMain(end+1)=it(4);
                TotalSolvTimeFromSolver(end+1)=it(5);
                DistDataTime(end+1)=it(6);
                testing(end+1)=it(7);
                timeMisc(end+1)=it(8);
                NLOPTTime(end+1)=it(9);
                CommuniTime(end+1)=it(10);
                scaledCommuniTime(end+1)=it(10)/it(1);
            end
            RegularDEgrees=         ['Degree:         ', num2str(G2)];
            disp(RegularDEgrees);
            IterRegular=            ['Iteration:      ',num2str(iterReg)];
            disp(IterRegular);
            TimeTotal=              ['Master   :      ',num2str(MasterTimeFromMain)];
            disp(TimeTotal);
            TotalTrain=             ['Train_Time:     ',num2str(TotalTrainTimeFromMain)];
            disp(TotalTrain);
            SolverTimeMain=         ['Solver_Main:    ',num2str(TotalSolvTimeFromMain)];
            disp(SolverTimeMain);
            SolverTimeSolver=       ['Solver_solver:  ',num2str(TotalSolvTimeFromSolver)];
            disp(SolverTimeSolver);
            DstrData=               ['Distr Data:     ',num2str(DistDataTime)];
            disp(DstrData);
            Test=                   ['Testing:        ',num2str(testing)];
            disp(Test);
            TimeMiscc=              ['Misc in main:   ',num2str(timeMisc)];
            disp(TimeMiscc)
            TimeNLOPT=              ['NLOPT time:     ',num2str(NLOPTTime)];
            disp(TimeNLOPT);
            TimeCommu=              ['Communi time:   ',num2str(CommuniTime)];
            disp(TimeCommu);
            TimeScaledCommunication=['Commui per Iter:' , num2str(scaledCommuniTime)];
            disp(TimeScaledCommunication);
        end
        %%%%%%%%%% Figures with all times
        [colorStyl]=legendNaedOnType(type, r, colors);
        
        figure(1)
        Plotting(figure(1),G2, iterReg, colorStyl, 'Iterations of Regular Graphs', legendOffFig, fontSize, 'Iterations','northeast')

        figure(2)
        Plotting(figure(2),G2, log(TotalTrainTimeFromMain), colorStyl, 'Training Time', legendOffFig, fontSize, 'log(Time)','northwest')

        figure(3)
        Plotting(figure(3),G2, scaledCommuniTime, colorStyl, 'Scaled Communication Time',  legendOffFig, fontSize, 'Time in sec','northwest')
        
    end
    r=r+1;
end
    SavingPlots(figure(1), 'IterationsRegVSMean')
    SavingPlots(figure(2), 'TrainingTime')
    SavingPlots(figure(3), 'ScaledCommuniTime')
    
function [ColorStyle]=legendNaedOnType(iType, counter,ColorVec)
if iType==1
    ColorStyle=['-',ColorVec{counter},'*'];
elseif iType==2
    ColorStyle=['--',ColorVec{counter+1},'o'];
elseif iType==3
    ColorStyle='-ro';
end
end
function SavingPlots(f, FileName)
set(f,'PaperOrientation','landscape');
set(f, 'PaperPosition',[0 0 11.5 8.5]);
set(gca,'FontSize',25);
ParentDir='.\..\pictures/';
fn=fullfile(ParentDir);
if ~exist(fn, 'dir')
    disp('Creating the directory for saving the plots ... !!')
    mkdir(fn)
end
Name1=[ParentDir,'\',FileName,'.eps'];
saveas(f, Name1,'epsc');
Name2=[ParentDir,'\',FileName,'.pdf'];
saveas(f, Name2);
%close(f)
end

function Plotting(fig,G2, FuncToPlot, colorStype, PlotTitle, legendOffFig, LegendFontSize, ylabelName,LegendLocation)
plot(G2,FuncToPlot,colorStype);
title(PlotTitle)
legend(legendOffFig(2:end),'Location', LegendLocation,'FontSize',LegendFontSize);
%legend(legendOffFig(2:end),'FontSize',LegendFontSize);
xlabel('Degree of graphs','FontSize',LegendFontSize-15);
ylabel(ylabelName,'FontSize',LegendFontSize);
xticks(G2);
xticklabels(G2);
set(gca,'FontSize',LegendFontSize);
% a = get(gca,'XTickLabel');
% set(a,'FontSize',LegendFontSize)
% ax=get(gca,'Xlabel');
% set(ax,'FontSize', Fontsize)
hold on
end

