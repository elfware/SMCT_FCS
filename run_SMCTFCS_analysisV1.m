%% run_SMCTFS_analysisV1
%This script analyses SMCT-FCS experimental data as in Marklund et al, 
%Nature, 2020. The script starts with rawdata from the SMCT-FCS microscope, 
%sort the datastreams, build trajectories, calculates auotcorelation
%functions for the individual trajectories.
%Then the script selects the sliding trajectores, pools the data for 
%different expermental days and, plots the mean of the autocorrelation 
%functions for individual trajectoreries. After calculating diffusion constanst for sliding,
%the script also fits a pitch dependent autocorerelation model
%to the data, which gives an estimate of the pitch for rotation coupled
%sliding
%
%Raw data is availible at elflab.icm.uu.se
%
%2020-04-05 Emil v1: First version of the streamlined pipeline for full 
%analysis. The script asumes that you have created a folder 'saved_files'
%in the same folder as you have the git repository,  as well that you have  raw data in a
%folder 'raw_data_SMCT_FCS', also in the same folder as were you have the
%git repository. You should thus have three folders, 'saved_files', 'raw_data_SMCT_FCS'
%and 'SMCT_FCS' from github in the same folder, which is one level above
%this script in the folder structure.
%On a single core the full analysis will take some hours to run
%(less than 24 hours)

%% Add paths
addpath(genpath('Functions/FunctionsConfocal'))
addpath(genpath('Functions/Polarization modeling'));

%Path to where you have the raw data. Raw data is availible at
%elflab.icm.uu.se
rawDataPath='../raw_data_SMCT_FCS';
s=what(rawDataPath);
rawDataPathAbs=s.path;
 
%path to were you are saving files. This folder must exist to be able to
%run the script.
savePath='../saved_files';
s=what(savePath);
savePathAbs=s.path;
%% Define folder names with data

%The names of the scripts in the analyseRawDataScripts folder, that analyses raw data for
%individual experiment days, should be on the form ['analyseRawData_'
%rawDataSets{i}]

%LacI-MBP-R
rawDataSets{1}='170705_LacIdim2MBP';
rawDataSets{2}='170713_LacIdim2MBP';
rawDataSets{3}='170718_LacIdim2MBP';
rawDataSets{4}='181010_LacIdim2MBP';
rawDataSets{5}='181021_LacIdim2MBP';
rawDataSets{6}='181105_LacIdim2MBP';


%LacI-R
rawDataSets{7}='170607_LacIdim2';
rawDataSets{8}='170704_LacIdim2';
rawDataSets{9}='180704_LacIdim2';
rawDataSets{10}='181105_LacIdim2';

%% Analyse rawdata
%Sort datastreams, build trajectories and calculate autocorrelation functions

for thisind=1:numel(rawDataSets)
    currDataPath=[rawDataPathAbs '/' rawDataSets{thisind}];
    run(['analyseRawDataScripts/analyseRawData_' rawDataSets{thisind}]);
    save([savePathAbs '/rawDataAnalysed_' rawDataSets{thisind}])
    clear Trajectories;
    clear Trajobject;
    clear Scanobject;
    clear Scanobjects;
end


close all;

%% %%
% Paths to analysed partially analysed data
pathAutoCorr=cell(1,numel(rawDataSets));
for thisind=1:numel(rawDataSets)
    pathAutoCorr{thisind}=[savePathAbs '/rawDataAnalysed_' rawDataSets{thisind}];
end

%Vector with indeces pointing to what rawdata path to use. Sliding and stationary molecules
%are picked out from the same analysed rawdatapaths
setToExps=[1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10];
numexperiments=numel(setToExps);
numSets=numel(pathAutoCorr);

%'Pos' = pick out sliding molecules, 'Stuck' = pick out stationary
%molecules
control={'Pos','Pos','Pos','Pos','Pos','Pos','Pos','Pos','Pos','Pos','Stuck','Stuck','Stuck','Stuck','Stuck','Stuck','Stuck','Stuck','Stuck','Stuck'};



normNumber=2; %Number of timelags to use in the normalization of the ACF.
integrateStart=1;

thelag=5;
stdCutOff=0.6;%Throw away molecule ACFs with a std larger than this. Effectively removes noisy and uninformative ACFs
stdCountCutOff=inf;%
maxMeanCount=inf;
minMeanCount=45;%Minimum mean count of the traces to be used for aucorrelation analysis. Also removes noisy autocorrelation function that would bias the mean otherwise.
minTotCounts=0;
%Cell structuring determining how data is pooled
%[1 2 3 4 5 6] is LacI-MBP-R sliding
%[7 8 9 10] is LacI-R sliding
%[11 12 13 14 15 16] is LacI-MBP-R stationary
%[17 18 19 20] is LacI-R stationary
combineExperiments={[1 2 3 4 5 6],[7 8 9 10],[11 12 13 14 15 16],[17 18 19 20]};
numCombines=numel(combineExperiments);



endInd=22;
timeWind=7;%Number of time lags in the autocorrelation function


shiftperiod=500e-6;%Laser on time
shifttime=166e-6;%laser off time. Dark time of the laser has been removed from autocorrelation function (first happens after 500 us)
%% Remove trajectories that are bleaching in more than one step
%manually picked out for the example dataset. here you load the indexes for those trajectories. If you change classification of
%sliding below, you cannot load those indexes! In that case you again have to pick out
%bleachers first! You can also set classparams.manualBleachRemoval to 0 and
%skip this step
classparams.manualBleachRemoval=1;

notBleaching=cell(1,numel(setToExps));
numNotBleaching=zeros(1,numel(setToExps));
numTotB=zeros(1,numel(setToExps));
if classparams.manualBleachRemoval==1
    load('notBleachingTrajs.mat','notBleaching');
end
    
%% Classification of sliding
clear diff;
load(pathAutoCorr{2},'Trajobject');
deltats=diff(Trajobject.TrajPosEstimation.FixWindowAv{1}(:,1));
deltat=deltats(1)*1000;%time step in seconds times 1000 to get in milliseconds
% Parameters for the sldiing classification. Find particles diffusing along
% the DNA axis
classparams.dE=deltat; %exposure time
classparams.nDeltaT=1;%steppeing size in deltat. if = 1 all time positions are used, if 2 everyother is skipped
classparams.deltat=classparams.nDeltaT*classparams.dE;
classparams.scale=1000; %nm/pixel, one pixel is 1 um
classparams.varwindow=10;%10
classparams.minPrincRelVar=0.75;%0.8 
classparams.minPrincAbsVar=0.0015;%1.5*0.05*deltat/200;%2
classparams.maxPrincAbsVar=inf;
classparams.minNotStuck=0.02*deltat/200;%0.5 for PCAver9
classparams.minMinWindowVar=0.00008;
classparams.maxMaxWindowVar=0.2*(deltat/200)*(classparams.varwindow/5);%
classparams.maxMinWindowVar=0.00035;
classparams.maxMeanDiff=inf;
% error estimation ffrom gain
classparams.minMeanDiff=-inf;
classparams.maxAbsVar=7e-4;
classparams.angledomain=[-pi/4 pi/4+pi/2];
minTrajTime=0.55;
classparams.startTime=0.2;
classparams.minTrajLength=round(minTrajTime*1000/deltat); %
classparams.msdStartPoint=1;
classparams.msdEndPoint=2;
classparams.minAngle=0-0.35;
classparams.maxAngle=0+0.35;
classparams.maxDiffusionsPrinc1=0.07;
classparams.minDiffusionsPrinc1=-0.005;
classparams.newDeltat=16;


diffMaxVec=0.075*ones(1,numel(setToExps));
diffMinVec=-inf*ones(1,numel(setToExps));

%diffMaxVec(combineExperiments{1})=0.02;
%diffMinVec(combineExperiments{1})=-0.005;

%diffMaxVec(combineExperiments{2})=0.07;
%diffMinVec(combineExperiments{2})=0.02;

%slidingcontrol='Neg';%IF Neg, slidingInd is set to alal trajectoreis outside of minAngle of maxAngle, that is the negative control is analysed for aoutocorelation; If stuck, pick out stuck trajectories
DataSelectionRules.Trajectory.Mode='FixWindowAv';%Use: FixWindowAv, RollingAv or DirectEstimation3Points
%sliding=cell(1,numexperiments);
%sliding=cell(1,numSets);
parfor theSet=1:numel(setToExps)
    classparamscurrent=classparams;
    classparamscurrent.maxDiffusionsPrinc1=diffMaxVec(theSet);
    classparamscurrent.minDiffusionsPrinc1=diffMinVec(theSet);
    disp(['sliding classidication for set ' num2str(theSet)])
    TrajobjectPar=load(pathAutoCorr{setToExps(theSet)},'Trajobject');
    [slidingInd,notstuck,angles,rellat,trajPrinclatent,thescores,okLength,trajL,classifiedParams] = slidingClassifyPCACameraTemplatever12_11(TrajobjectPar.Trajobject,DataSelectionRules,classparamscurrent,1);
    allInds=1:numel(notstuck);
    classparams2=classparamscurrent;

    classparams2.maxAngle=+inf;
    classparams2.minAngle=-inf;
    [slidingIndalldir,~,~,~,~,~] = slidingClassifyPCACameraTemplatever12_11(TrajobjectPar.Trajobject,DataSelectionRules,classparams2,1);
    
    if strcmp(control{theSet},'Neg')
        slidingInd=setdiff(slidingIndalldir,slidingInd);
    end
    
    if strcmp(control{theSet},'Stuck')
        slidingInd=allInds(~notstuck);
 %       slidingInd=slidingInd(trajL(slidingInd)>classparams.minTrajLength);
    end
    
    
    
    
    disp(['For set ' num2str(theSet) ', Number of trajectories longer than treshhold: ' num2str(okLength)]);
    disp(['For set ' num2str(theSet) ', Number of trajectories passing classification: ' num2str(numel(slidingInd))]);
%       if ~strcmp(trajIndexVarName{theSet},'slidingInd')
%         obj=MMload(pathTrajIndex{theSet},trajIndexVarName{theSet});
%         slidingInd=obj.(trajIndexVarName{theSet});
%           
%       end
%% Remove trajectories that bleaches in in two steps or more, here done manually, remove this section if sliding classification params are changed!!
slidingIndPreBleachRem=slidingInd;
if classparams.manualBleachRemoval
    if ~isempty(notBleaching{theSet})
        numTotB(theSet)=numel(slidingInd);
        slidingInd=slidingInd(notBleaching{theSet});
        numNotBleaching=numel(slidingInd);
    end
end
    %%
    slidingClassified=classifiedParams(slidingInd,:);
    
    %pick out sliding trajectories to its own variable
    slidingPixelScores=cell(1,numel(slidingInd));
    slidingTrajs=cell(1,numel(slidingInd));
    startInd=ceil(classparams.startTime/(deltat*0.001));
    trajTime=zeros(1,numel(slidingInd));
    for i=1:numel(slidingInd)
        slidingPixelScores{i}=thescores{slidingInd(i)}(1:classparams.nDeltaT:end,:);
        slidingTrajs{i}=TrajobjectPar.Trajobject.TrajPosEstimation.FixWindowAv{slidingInd(i)}(startInd:classparams.nDeltaT:end,2:3);
        trajTime(i)=TrajobjectPar.Trajobject.Traj(slidingInd(i),1); %Total time in s
    end
    
    
    %scale=1000; %tthescore and slidingPixelsScors is in units of um, so 1 pixel is 1 um, scale is in of nm/pixel, so scale=1000
    Dimension=2;
    [avgTotDiff,compTotDiff]=getTotalDiffCVE(slidingPixelScores,classparams.deltat,classparams.scale,Dimension);
    [avgDiff,compDiff]=getDiffCVECamera(slidingPixelScores,classparams.deltat,classparams.scale,Dimension);
    
    minUseLength=round((minTrajTime-classparams.startTime)*1000/deltat)-1;
    
   [msd,compDiffMSD]=getDiffMSDCamera_fixdiv2_ver2(slidingPixelScores,classparams.deltat,classparams.scale,Dimension,classparams.msdStartPoint,classparams.msdEndPoint);
   Dxs=compDiff(1,:);
   Dys=compDiff(2,:);
    
    DxsMSD=compDiffMSD(1,:);
    DysMSD=compDiffMSD(2,:);
        
    
    varPrinc=trajPrinclatent(1,slidingInd);

%     
     disp(['For set ' num2str(theSet) ', Number of trajectories slidng in all directions: ' num2str(numel(slidingIndalldir))]);
     
     angleTot=classparams.maxAngle-classparams.minAngle;
     falsePos=((numel(slidingIndalldir)-numel(slidingIndPreBleachRem))/(pi-angleTot))/(numel(slidingIndPreBleachRem)/angleTot);
     disp(['For set ' num2str(theSet) 'Estimated fraction false positive sliders :' num2str(falsePos)]);
     
     [diffPrinc1,diffPrinc2,allDiffPrinc1,allDiffPric2,meanDiffPrinc1,meanDIffPrinc2]=calcTrajStepdiff(thescores,1,1,slidingInd);
     
    
    
     
     sliding{theSet}.numTotTrajs=numel(notstuck);
     sliding{theSet}.notsutck=notstuck;
     sliding{theSet}.stuckInd=allInds(~notstuck);
     sliding{theSet}.numTrajsOkLength=okLength;
     sliding{theSet}.slidingInd=slidingInd;
     sliding{theSet}.slidingIndalldir=slidingIndalldir;
     sliding{theSet}.falsePos=falsePos;
     sliding{theSet}.Dxs=Dxs;
     sliding{theSet}.Dys=Dys;
     sliding{theSet}.DxsMSD=DxsMSD;
     sliding{theSet}.DysMSD=DysMSD;
     
     sliding{theSet}.varPrinc=varPrinc;
     sliding{theSet}.meanDiffPrinc1=meanDiffPrinc1(slidingInd);
     sliding{theSet}.slidingTrajs=slidingTrajs;
     sliding{theSet}.slidingPixelScores=slidingPixelScores;
     sliding{theSet}.msd=msd;
     sliding{theSet}.trajTime=trajTime; %total time in s of slidign trajectories, before cutoff first part of trajecotry
     sliding{theSet}.slidingClassified=slidingClassified;
     sliding{theSet}.slidingIndPreBleachRem=slidingIndPreBleachRem;
    
end



%% Calcualte diffusion constants
clear diff;
newDeltat=16;%classparams.newDeltat;
kWind=floor(newDeltat/deltat);
newStartPoint=1;
newEndPoint=10;
slidingAvg=cell(1,30);
%newmaxDiffusionsPrinc1=inf;
for theSet=1:numel(setToExps)
     disp(['Calculating diffusion for set ' num2str(theSet)])
    slidingPixelScores=sliding{theSet}.slidingPixelScores;
    numtrajs=numel(slidingPixelScores);
    averagePos=cell(1,numtrajs);
    allDNADisplacements=[];
    for i=1:numtrajs
        if ~isempty(slidingPixelScores{i})
        smoothAvg1=smooth(slidingPixelScores{i}(:,1),kWind);
        smoothAvg2=smooth(slidingPixelScores{i}(:,2),kWind);
        averagePos{i}=[];
        averagePos{i}(:,1)=smoothAvg1(ceil(kWind/2):kWind:end);
        averagePos{i}(:,2)=smoothAvg2(ceil(kWind/2):kWind:end);
        %averagePos{i}(:,1)=smoothAvg1(1:kWind:end);
        %averagePos{i}(:,2)=smoothAvg2(1:kWind:end);
        displacements=diff(averagePos{i}(:,1));
        allDNADisplacements=[allDNADisplacements;displacements];
        end
    end
    
    Dimension=2;
    [avgTotDiff,compTotDiff]=getTotalDiffCVE(averagePos,newDeltat,classparams.scale,Dimension);
    [avgDiff,compDiff]=getDiffCVECamera(averagePos,newDeltat,classparams.scale,Dimension);
    [msd,compDiffMSD]=getDiffMSDCamera_fixdiv2_ver2(averagePos,newDeltat,classparams.scale,Dimension,newStartPoint,newEndPoint);
        
    Dxs=compDiff(1,:);
    Dys=compDiff(2,:);
    
    DxsMSD=compDiffMSD(1,:);
    DysMSD=compDiffMSD(2,:);
  
    %skew=skewness(allDNADisplacements);
    n=numel(allDNADisplacements);
    errorskew=sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3)));
    
    %disp(['skewness : ' num2str(skew) ' +/- ' num2str(errorskew) ' um']);
    
    slidingAvg{theSet}.Dxs=Dxs;
    slidingAvg{theSet}.Dys=Dys;
    slidingAvg{theSet}.DxsMSD=DxsMSD;
    slidingAvg{theSet}.DxsMSD=DxsMSD;
    slidingAvg{theSet}.slidingInd=sliding{theSet}.slidingInd;
    slidingAvg{theSet}.compTotDiff=compTotDiff;
    slidingAvg{theSet}.averagePos=averagePos;
    slidingAvg{theSet}.allDNADisplacements = allDNADisplacements;
    
   
    
    
    
end
%% PLot MSD plots for sliders
%figure();hist(compTotDiffboot(1,:));
theStartP=1;
theEndP=3;
averagePoses1=[slidingAvg{1}.averagePos slidingAvg{2}.averagePos slidingAvg{3}.averagePos slidingAvg{4}.averagePos slidingAvg{5}.averagePos slidingAvg{6}.averagePos];
averagePoses2=[slidingAvg{7}.averagePos slidingAvg{8}.averagePos slidingAvg{9}.averagePos slidingAvg{10}.averagePos];
[msd2Tot,compDiffMSD2Tot,msdStd2Tot,msdN2Tot]=getTotDiffMSDCamera_fixdiv2_ver3(averagePoses2,newDeltat,classparams.scale,Dimension,theStartP,theEndP);
[msd1Tot,compDiffMSD1Tot,msdStd1Tot,msdN1Tot]=getTotDiffMSDCamera_fixdiv2_ver3(averagePoses1,newDeltat,classparams.scale,Dimension,theStartP,theEndP);

[msd2,compDiffMSD2]=getDiffMSDCamera_fixdiv2_ver2(averagePoses2,newDeltat,classparams.scale,Dimension,theStartP,theEndP);
[msd1,compDiffMSD1]=getDiffMSDCamera_fixdiv2_ver2(averagePoses1,newDeltat,classparams.scale,Dimension,theStartP,theEndP);


tvals=newDeltat:newDeltat:newDeltat*(theEndP-theStartP+1);


figure();
errorbar(tvals,msd1Tot{1},msdStd1Tot{1}./sqrt(msdN1Tot{1}),'bd')
hold on;
errorbar(tvals,msd2Tot{1},msdStd2Tot{1}./sqrt(msdN2Tot{1}),'rd')
[P1,S1]=polyfit(tvals,msd1Tot{1},1);
plot([tvals(1)*0.5 tvals tvals(end)*1.2],polyval(P1,[tvals(1)*0.5 tvals tvals(end)*1.2]),'b');

[P2,S2]=polyfit(tvals,msd2Tot{1},1);
plot([tvals(1)*0.5 tvals tvals(end)*1.2],polyval(P2,[tvals(1)*0.5 tvals tvals(end)*1.2]),'r');
xlabel('Time (ms)');
ylabel('MSD micronsquared');
%

figure();
errorbar(tvals,mean(msd1{1}),std(msd1{1})./sqrt(numel(msd1{1})),'bd')
hold on;
errorbar(tvals,mean(msd2{1}),std(msd2{1})./sqrt(numel(msd2{1})),'rd')
[P1,S1]=polyfit(tvals,mean(msd1{1}),1);
plot([tvals(1)*0.5 tvals tvals(end)*1.2],polyval(P1,[tvals(1)*0.5 tvals tvals(end)*1.2]),'b');

[P2,S2]=polyfit(tvals,mean(msd2{1}),1);
plot([tvals(1)*0.5 tvals tvals(end)*1.2],polyval(P2,[tvals(1)*0.5 tvals tvals(end)*1.2]),'r');
xlabel('Time (ms)');
ylabel('MSD micronsquared ');
%%


numTimeScales=endInd-timeWind+1;
scaleRes=cell(1,numTimeScales);
timeRegimes=zeros(1,numTimeScales);
combCorrTimesMean=zeros(numel(combineExperiments),numTimeScales);
resCorrTimesMean=zeros(numel(setToExps),numTimeScales);
combCorrTimesStd=zeros(numel(combineExperiments),numTimeScales);
resCorrTimesStd=zeros(numel(setToExps),numTimeScales);

combCorrTimesLogMean=zeros(numel(combineExperiments),numTimeScales);
resCorrTimesLogMean=zeros(numel(setToExps),numTimeScales);
combCorrTimesLogStd=zeros(numel(combineExperiments),numTimeScales);
resCorrTimesLogStd=zeros(numel(setToExps),numTimeScales);

resCorrTimesError=zeros(numel(setToExps),numTimeScales);

combCorrTimesMeanUpper=zeros(numel(combineExperiments),numTimeScales);
combCorrTimesMeanLower=zeros(numel(combineExperiments),numTimeScales);

combCorrTimesStdUpper=zeros(numel(combineExperiments),numTimeScales);
combCorrTimesStdLower=zeros(numel(combineExperiments),numTimeScales);

combCorrTimesError=zeros(numel(combineExperiments),numTimeScales);


corrTimesMatr=cell(1,numel(setToExps));

resCorrTimesAll=cell(1,numel(setToExps));
T=[];

%%
for theExp=cell2mat(combineExperiments)
    disp(['Loading data for experiment ' num2str(theExp)]);
    Trajobjects{theExp}=load(pathAutoCorr{setToExps(theExp)},'Trajobject');
end
%% Calculate and filter out individiuval trajectory auto correlationfunctions
%Also pool data into the combines structures, to compare all data for each protein 
for timescaleInd=4:17%4:numTimeScales
    minLagInd=timescaleInd;
    %maxLagInd=endInd;
    
    endNormInd=minLagInd+timeWind-1;
    maxLagInd=endNormInd;
   
    disp(['Time ' num2str(timescaleInd) ' out of ' num2str(numTimeScales)])
    if ~isempty(T)
            disp(['Start of tiem interval: ' num2str(T(timescaleInd))]);
    end
    %%
    for theExp=cell2mat(combineExperiments)%[combineExperiments{1} combineExperiments{2} combineExperiments{3} combineExperiments{4} combineExperiments{5} combineExperiments{6} combineExperiments{7} combineExperiments{8} combineExperiments{9} combineExperiments{10} combineExperiments{11}]
        disp(['Analyzing experiment ' num2str(theExp)])
        
        Trajobject=Trajobjects{theExp}.Trajobject;
        
        %data=load(pathTrajIndex{theExp},trajIndexVarName{theExp});
        %nameIndexVar=fieldnames(data);
        trajIndexOfInterest=sliding{theExp}.slidingInd;
        numOfInterest=length(trajIndexOfInterest);
        
        RChSum=Trajobject.TTTTRTrajs.Autocorrelation.RChSum;
        RChCross=Trajobject.TTTTRTrajs.Autocorrelation.RChCross;
        Timechift=Trajobject.TTTTRTrajs.Autocorrelation.Timechift;
        
        
        T=mean(Timechift(trajIndexOfInterest,:),1);
        % Fix times in autocorrelation to  correct for the shift
        T=T+floor(T/shiftperiod)*shifttime;
        maxLagTime=T(maxLagInd);
        plotT=0:1e-7:maxLagTime;
        
        %         minLagInd=find(T>minLagTime,1);
        %         maxLagInd=find(T>maxLagTime,1)-1;
        %         endNormInd=find(T>endNormLagTime,1)-1;
        CorrSum=sum(RChSum(trajIndexOfInterest,:));
        
        
        autoCorrValsAllRaw=RChSum(trajIndexOfInterest,minLagInd:maxLagInd);
        autoCrossCorrValsAllRaw=RChCross(trajIndexOfInterest,:);
        
        autoCorrValsAll=autoCorrValsAllRaw;
        initialAmplitude=autoCorrValsAll(:,1);
        autoCorrValsAll=autoCorrValsAll-autoCorrValsAll(:,maxLagInd-minLagInd+1)*ones(1,numel(minLagInd:maxLagInd));
        autoCorrValsAll=autoCorrValsAll./abs(mean(autoCorrValsAll(:,1:normNumber),2)*ones(1,numel(minLagInd:maxLagInd)));
        
    
        newNum=numel(autoCorrValsAllRaw(1,thelag+1:end));
        endPoints=autoCorrValsAllRaw(:,end)*ones(1,newNum);
        
        midPoints=autoCorrValsAllRaw(:,thelag+1:end)-endPoints;
        startPoints=autoCorrValsAllRaw(:,1:end-thelag)-endPoints;
        thediffs=midPoints-startPoints;
   

        autoCorrValsAllQ=midPoints./startPoints;
        
        autoCorrValsAllQ(midPoints<0)=0;
        autoCorrValsAllQ(midPoints>startPoints)=1;
        autoCorrValsAllQ(startPoints<0)=NaN;

        
        
        autoCorrStd=std(autoCorrValsAll(:,:),0,2);
        
        
        totCountsAll=sum(Trajobject.Traj(trajIndexOfInterest,6:7),2);
        meanPolAll=Trajobject.Traj(trajIndexOfInterest,4);
        countsAll=cell(1,numOfInterest);
        meanCountsAll=zeros(numOfInterest,1);
        stdCountsAll=zeros(numOfInterest,1);
        k=1;
        for theIndex=trajIndexOfInterest
            countsAll{k}=sum(Trajobject.TrajPosEstimation.FixWindowAv{theIndex}(:,6:7),2);
            stdCountsAll(k)=std(countsAll{k});
            meanCountsAll(k)=mean(countsAll{k});
            k=k+1;
        end
        
        okAutoCorrs=autoCorrStd<stdCutOff;
        okAutoCorrs=(okAutoCorrs+((stdCountsAll./meanCountsAll)<stdCountCutOff))>1;
        okAutoCorrs=(okAutoCorrs+(autoCorrValsAll(:,1)>0))>1;
        okAutoCorrs=(okAutoCorrs+(totCountsAll>minTotCounts))>1;
        
        okAutoCorrs=(okAutoCorrs+(meanCountsAll<maxMeanCount))>1;
        okAutoCorrs=(okAutoCorrs+(meanCountsAll>minMeanCount))>1;
        %okAutoCorrs=(okAutoCorrs+(autoCrossCorrValsAllRaw(:,1)==0))>1;
        
        numCrossSingle=sum(autoCrossCorrValsAllRaw(:,1)==0);
        % okAutoCorrs=results{theExp}.okAutoCorrs;
        %okAutoCorrs=chooseInd.results{theExp}.okAutoCorrs;
        
        slidingIndok=trajIndexOfInterest(okAutoCorrs);
        numOkAutoCorrs=sum(okAutoCorrs);
        autoCorrValsAllok=autoCorrValsAll(okAutoCorrs,:);
        autoCorrValsAllQok=autoCorrValsAllQ(okAutoCorrs,:);
        autoCorrValsAllRawok=autoCorrValsAllRaw(okAutoCorrs,:);
        autoCorrVals=mean(autoCorrValsAllok,1);
        lagTimes=T(minLagInd:maxLagInd);

        
        newEndNormInd=endNormInd-minLagInd+1;
        
 
        
        corrTimesAll=cumtrapz(lagTimes(integrateStart:newEndNormInd),autoCorrValsAll(:,integrateStart:newEndNormInd)');
        corrTimesAllok=cumtrapz(lagTimes(integrateStart:newEndNormInd),autoCorrValsAllok(:,integrateStart:newEndNormInd)');
        
        
        
        
        
        
        
        results{theExp}.initialAmplitude=initialAmplitude;
        results{theExp}.lagTimes=lagTimes;
        results{theExp}.autoCorrValsAll=autoCorrValsAll;
        results{theExp}.okAutoCorrs=okAutoCorrs;
        results{theExp}.autoCorrValsAllok = autoCorrValsAllok;
        results{theExp}.corrTimesAllok=corrTimesAllok;
        results{theExp}.corrTimesAll=corrTimesAll;
        results{theExp}.totCountsAll=totCountsAll;
        results{theExp}.meanPolAll=meanPolAll;
        results{theExp}.slidingIndok=slidingIndok;
        results{theExp}.autoCorrStd=autoCorrStd;
        results{theExp}.autoCorrValsAllQ=autoCorrValsAllQ;
        results{theExp}.autoCorrValsAllQok=autoCorrValsAllQok;
        results{theExp}.autoCorrValsAllRaw= autoCorrValsAllRaw;
        results{theExp}.countsAll=countsAll;
        results{theExp}.stdCountsAll=stdCountsAll;
        results{theExp}.meanCountsAll=meanCountsAll;
        results{theExp}.numCrossSingle=numCrossSingle;
        
        
    end
    
    %% Combine/ pool data from same type of experiment together
    %combines=cell(1,numCombines);
    %cbootStrap=cell(1,numCombines);
    for theComb=[1:numel(combineExperiments)]
        disp(['Pooling data for combination number ' num2str(theComb)])
%         cols=numel(results{combineExperiments{theComb}(1)}.autoCorrVals);
%         autoCorrNormVals=zeros(1,cols);
        amplitudes=ones(1,numel(combineExperiments{theComb}));
        corrTimesAll=[];
        corrTimesAllok=[];
        varPrincsok=[];
        totCountsAllok=[];
        totCountsAll=[];
        Dxs=[];
        Dxsok=[];
        DxsMSDok=[];
        DxsMSD=[];
        meanPolAllok=[];
        meanPolAll=[];
        autoCorrValsAllok=[];
        meanDiffPrinc1=[];
        autoCorrStd=[];
        autoCorrValsAllQ=[];
        autoCorrValsAllQok=[];
        autoCorrValsAllRaw=[];
        for theSet=1:numel(combineExperiments{theComb})
            theExp=combineExperiments{theComb}(theSet);
            %amplitudes(theSet)=results{theExp}.initialAmplitude*(results{theExp}.fitNorm-results{theExp}.endNorm);
            corrTimesAllok=[corrTimesAllok results{theExp}.corrTimesAllok];
            corrTimesAll=[corrTimesAll results{theExp}.corrTimesAll];
            varPrincsok=[varPrincsok sliding{theExp}.varPrinc(results{theExp}.okAutoCorrs)];
            totCountsAllok=[totCountsAllok;results{theExp}.totCountsAll(results{theExp}.okAutoCorrs)];
            totCountsAll=[totCountsAll;results{theExp}.totCountsAll];
            Dxs=[Dxs sliding{theExp}.DxsMSD];
            Dxsok=[Dxsok sliding{theExp}.DxsMSD(results{theExp}.okAutoCorrs)];
            DxsMSDok=[DxsMSDok sliding{theExp}.DxsMSD(results{theExp}.okAutoCorrs)];
            DxsMSD=[DxsMSD sliding{theExp}.DxsMSD];
            meanPolAllok=[meanPolAllok;results{theExp}.meanPolAll(results{theExp}.okAutoCorrs)];
            meanPolAll=[meanPolAll;results{theExp}.meanPolAll];
            autoCorrValsAllok=[autoCorrValsAllok;results{theExp}.autoCorrValsAllok];
            meanDiffPrinc1=[meanDiffPrinc1 sliding{theExp}.meanDiffPrinc1];
            autoCorrStd=[autoCorrStd;results{theExp}.autoCorrStd];
            autoCorrValsAllQ=[autoCorrValsAllQ;results{theExp}.autoCorrValsAllQ];
            autoCorrValsAllQok=[autoCorrValsAllQok;results{theExp}.autoCorrValsAllQok];
            autoCorrValsAllRaw=[autoCorrValsAllRaw;results{theExp}.autoCorrValsAllRaw];
            
        end
        medianDxsMSDok=median(DxsMSDok);
    
        
        
      
        combines{theComb}.corrTimesAll= corrTimesAll;
        combines{theComb}.varPrincsok =varPrincsok;
        combines{theComb}.Dxsok=Dxsok;
        combines{theComb}.DxsMSDok=DxsMSDok;
        combines{theComb}.DxsMSD=DxsMSD;
        combines{theComb}.meanPolAllok=meanPolAllok;
        combines{theComb}.meanPolAll=meanPolAll;
        combines{theComb}.autoCorrValsAllok=autoCorrValsAllok;
        combines{theComb}.medianDxsMSDok=medianDxsMSDok;
        combines{theComb}.corrTimesAllok=corrTimesAllok;
  
        combines{theComb}.meanDiffPrinc1=meanDiffPrinc1;
        combines{theComb}.autoCorrStd=autoCorrStd;
     
        combines{theComb}.totCountsAll=totCountsAll;
        combines{theComb}.totCountsAllok = totCountsAllok;
        combines{theComb}.autoCorrValsAllQ=autoCorrValsAllQ;
        combines{theComb}.autoCorrValsAllQok=autoCorrValsAllQok;
        combines{theComb}.autoCorrValsAllRaw = autoCorrValsAllRaw;
        combines{theComb}.Dxs=Dxs;
        
    end
  
        
 

    numCombs=numel(combineExperiments);

    
    timeRegimes(timescaleInd)=lagTimes(1);
    scaleRes{timescaleInd}.lagTimes=lagTimes;
    
    scaleRes{timescaleInd}.combcorrTimesAllok=cell(1,numCombs);
    maxCorrTime(timescaleInd)=lagTimes(timeWind)-lagTimes(1);

    
    for i=1:numCombs
    
        scaleRes{timescaleInd}.combines{i}=combines{i};
        combCorrTimesMean(i,timescaleInd)=nanmean(combines{i}.corrTimesAllok(end,:));
        combCorrTimesStd(i,timescaleInd)=nanstd(combines{i}.corrTimesAllok(end,:));
        combCorrTimesError(i,timescaleInd)= combCorrTimesStd(i,timescaleInd)/sqrt(numel(combines{i}.corrTimesAllok(end,:)));
        inds=combines{i}.corrTimesAllok(end,:)>0;
        
        combCorrTimesLogMean(i,timescaleInd)=nanmean(log10(combines{i}.corrTimesAllok(end,inds)));
        combCorrTimesLogStd(i,timescaleInd)=nanstd(log10(combines{i}.corrTimesAllok(end,inds)));

        
        numRes=numel(combineExperiments{i});
        
         scaleRes{timescaleInd}.rescorrTimesAllok=cell(1,numRes);
        
        for j=1:numRes
            
            ind=combineExperiments{i}(j);
           % corrTimesMatr{ind}(:,timescaleInd)=results{ind}.corrTimesAllok(end,:)';
            scaleRes{timescaleInd}.results{ind}=results{ind};
           % scaleRes{timescaleInd}.results{ind}.okAutoCorrs=results{ind}.okAutoCorrs;
           % scaleRes{timescaleInd}.results{ind}.autoCorrValsAllok = results{ind}.autoCorrValsAllok;
           
            [~,numtrajectories]=size(results{ind}.corrTimesAll);
            if isempty(resCorrTimesAll{ind})
                resCorrTimesAll{ind}=zeros(numtrajectories,numTimeScales);
            end
                
            for theTraj=1:numtrajectories
                resCorrTimesAll{ind}(theTraj,timescaleInd)=results{ind}.corrTimesAll(end,theTraj);
            end
            
            resCorrTimesMean(ind,timescaleInd)=nanmean(results{ind}.corrTimesAllok(end,:));
            resCorrTimesStd(ind,timescaleInd)=nanstd(results{ind}.corrTimesAllok(end,:));
            resCorrTimesError(ind,timescaleInd)=resCorrTimesStd(ind,timescaleInd)./(sqrt(numel(results{ind}.corrTimesAllok(end,:))));
            
            inds=results{ind}.corrTimesAllok(end,:)>0;
              
            resCorrTimesLogMean(ind,timescaleInd)=nanmean(log10(results{ind}.corrTimesAllok(end,inds)));
            resCorrTimesLogStd(ind,timescaleInd)=nanstd(log10(results{ind}.corrTimesAllok(end,inds)));
            
        end

        
    end
    
 

end
clear Trajobjects;

%% Estimate rate/frequency of sliders appearing
clear diff;
for theSet=1:numel(setToExps);
    disp(['Event rate calculation for set' num2str(theSet)])
    TrajobjectPar=load(pathAutoCorr{setToExps(theSet)},'Trajobject');
    thediffs=diff(TrajobjectPar.Trajobject.Traj(:,8));
    inds=1:numel(TrajobjectPar.Trajobject.Traj(:,8));
    posStartInds=inds([1;thediffs]>0.5);
    measTime=0;
    for j=1:numel(posStartInds)-1
        measTime=measTime+TrajobjectPar.Trajobject.TrajPosEstimation.FixWindowAv{posStartInds(j+1)-1}(end,1)-TrajobjectPar.Trajobject.TrajPosEstimation.FixWindowAv{posStartInds(j)}(1,1);
        
    end
    
    measTime=measTime+TrajobjectPar.Trajobject.TrajPosEstimation.FixWindowAv{end}(end,1)-TrajobjectPar.Trajobject.TrajPosEstimation.FixWindowAv{posStartInds(end)}(1,1);
    sliding{theSet}.measTime=measTime;
end

slidingCombines=cell(1,numel(combineExperiments));

  for theComb=[1:numel(combineExperiments)]
        disp(['combination number ' num2str(theComb)])
        measTimeTot=0;
        numTrajsokLengthTot=0;
        numTotTrajsTot=0;
        numSlidingTot=0;
        for theSet=1:numel(combineExperiments{theComb})
            ind=combineExperiments{theComb}(theSet);
            measTimeTot=measTimeTot+sliding{ind}.measTime;
            numTrajsokLengthTot=numTrajsokLengthTot+sliding{ind}.numTrajsOkLength;
            numTotTrajsTot=numTotTrajsTot+sliding{ind}.numTotTrajs;
            numSlidingTot=numSlidingTot+numel(sliding{ind}.slidingInd);
        end
     slidingCombines{theComb}.measTime=measTimeTot;   
     slidingCombines{theComb}.numSliding=numSlidingTot;
     slidingCombines{theComb}.numTrajsOkLength=numTrajsokLengthTot;
     slidingCombines{theComb}.numTotTrajs=numTotTrajsTot;
  end      

  %%
  %% Plot autocorrealtion functions


times=12;%Time index of the first time lag of the autocorrelation you want to use

numplots=numel(times);
figure();
ind1=1;
ind2=2;
k=1;
perc=0;
numBoots=2000;
for timei=times
    
    [numtrajs1,~]=size(scaleRes{timei}.combines{ind1}.autoCorrValsAllok);
    [numtrajs2,~]=size(scaleRes{timei}.combines{ind2}.autoCorrValsAllok);
    bootCorrs1=zeros(numBoots,timeWind);
    bootCorrs2=zeros(numBoots,timeWind);
    
    for boot=1:numBoots
        bootInds1=ceil(rand(1,numtrajs1)*numtrajs1);
        bootInds2=ceil(rand(1,numtrajs2)*numtrajs2);
        
        bootCorrs1(boot,:)=trimmean(scaleRes{timei}.combines{ind1}.autoCorrValsAllok(bootInds1,:),perc);
        bootCorrs2(boot,:)=trimmean(scaleRes{timei}.combines{ind2}.autoCorrValsAllok(bootInds2,:),perc);
    
    
    end
    
    
    %subplot(numplots,1,k); 
    hold on
    errorbar(scaleRes{timei}.lagTimes,mean(bootCorrs1),std(bootCorrs1),'b');
    hold on;
    errorbar(scaleRes{timei}.lagTimes,mean(bootCorrs2),std(bootCorrs2),'r');
    k=k+1;
    %set(gca,'YScale','log');
    set(gca,'XScale','log');
    xlim([T(times(1)) T(times(end)+timeWind-1)])
    ylim([0 1.2]);
end
title('Sliding molecules');
xlabel('Time lag [s]');
ylabel('Autocorrelation');

numplots=numel(times);
figure();
ind1=3;
ind2=4;
k=1;
numBoots=2000;
for timei=times
    
    [numtrajs1,~]=size(scaleRes{timei}.combines{ind1}.autoCorrValsAllok);
    [numtrajs2,~]=size(scaleRes{timei}.combines{ind2}.autoCorrValsAllok);
    bootCorrs1=zeros(numBoots,timeWind);
    bootCorrs2=zeros(numBoots,timeWind);
    
    for boot=1:numBoots;
        bootInds1=ceil(rand(1,numtrajs1)*numtrajs1);
        bootInds2=ceil(rand(1,numtrajs2)*numtrajs2);
        
        bootCorrs1(boot,:)=trimmean(scaleRes{timei}.combines{ind1}.autoCorrValsAllok(bootInds1,:),perc);
        bootCorrs2(boot,:)=trimmean(scaleRes{timei}.combines{ind2}.autoCorrValsAllok(bootInds2,:),perc);
    
    
    end
    
    
    %subplot(numplots,1,k); 
    hold on
    errorbar(scaleRes{timei}.lagTimes,mean(bootCorrs1),std(bootCorrs1),'b');
    hold on;
    errorbar(scaleRes{timei}.lagTimes,mean(bootCorrs2),std(bootCorrs2),'r');
    k=k+1;
    %set(gca,'YScale','log');
    set(gca,'XScale','log');
    xlim([T(times(1)) T(times(end)+timeWind-1)])
    ylim([0 1.2]);
end
title('Stationary molecules')
xlabel('Time lag [s]');
ylabel('Autocorrelation');


%%

DxsMBP1=[slidingAvg{1}.Dxs slidingAvg{2}.Dxs slidingAvg{3}.Dxs];
DxsRhoda1=[slidingAvg{7}.Dxs slidingAvg{8}.Dxs];

DxsMBP2=[slidingAvg{4}.Dxs slidingAvg{5}.Dxs slidingAvg{6}.Dxs];
DxsRhoda2=[slidingAvg{9}.Dxs slidingAvg{10}.Dxs];



allDxs{1}=DxsMBP1;
allDxs{2}=DxsRhoda1;
allDxs{3}=DxsMBP2;
allDxs{4}=DxsRhoda2;

allDxs{1}=[DxsMBP1 DxsMBP2];
allDxs{2}=[DxsRhoda1 DxsRhoda2];


maxD=inf;
minD=-inf;

trims=0:80;
meansDxs=zeros(4,numel(trims));
errorsDxs=zeros(4,numel(trims));
meanBoots=cell(1,numel(trims));
numBoots=500;
meanBoot=zeros(numBoots,1);
k=1;
for i=trims
    meanBoots{k}=zeros(numBoots,size(meansDxs,1));
    disp(i);
    for j=1:numel(allDxs)
        currVals=allDxs{j}(allDxs{j}<maxD);
        currVals=currVals(currVals>minD);
        meansDxs(j,k)=trimmean(currVals,i);
        meansBoot=zeros(1,numBoots);
        numSamples=numel(currVals);
        
        for bootInd=1:numBoots
            sampleInds=ceil(rand(1,numSamples)*numSamples);
            meanBoot(bootInd)=trimmean(currVals(sampleInds),i);
        end    
        errorsDxs(j,k)=std(meanBoot);
        meanBoots{k}(:,j)=meanBoot;
    end
    k=k+1;
end

confInterval=0.95;
sortedFactors=sort(meanBoots{1}(:,2)./meanBoots{1}(:,1),'ascend');
minInd=round(numBoots*(1-confInterval)/2);
lowerBoundkFactor=sortedFactors(minInd);
upperBoundkFactor=sortedFactors(end-minInd);



figure();errorbar(trims,meansDxs(1,:),errorsDxs(1,:),'b');
hold on;
errorbar(trims,meansDxs(2,:),errorsDxs(2,:),'r');
ylabel('Estimated diffusion constant [um^2/s]');
xlabel('Trim percent in trimmean');
legend('LacI-MBP-R','LacI-R')
%% Fit autocorrelation model to the difference in auto correlation of the two proteins
totNumParams=2;
kguessExpMax=8;
kguessExpMin=1;

aguessExpMax=-6;
aguessExpMin=-7;

numGuesses=200;%Number of itnial guesses to to pe each fit (on data and bootstrapped data) 
numTerms=1;

%Index of timescales to use
%theScales=[5:26];
theScales=[timei];

numScales=numel(theScales);

%Datasets/combinations to compare 
compare{1}=[1,2];
bgInds=[3 4];

%Background or amplitude of background isn't taken into account in this version of the fit
combUnNorm{bgInds(1)}.autoCorrScaled=zeros(1,timeWind);
combUnNorm{bgInds(2)}.autoCorrScaled=zeros(1,timeWind);
combUnNorm{bgInds(1)}.startScale=timei;
combUnNorm{bgInds(2)}.startScale=timei;

D2=meansDxs(2,21); %Diffusion constants are here the 20% trimmed means
D=meansDxs(1,21);
kFactor=D2/D;

%set lower and Upper bounds for the parameters of the fit. In this version,
%the exponential k of LacI-MBP-R is te only true free parameter, and the
%the exponential for the other prother can vary within what is predicted
%from the 95% confidence itnerval of the translational diffusion constants
%of the protein
lb=zeros(1,3);
ub=zeros(1,3);
lb(1:totNumParams/2)=lowerBoundkFactor;%kFactor, so that k*kFactor is the exponantial for the second protein
lb(totNumParams/2+1:totNumParams)=-inf;% %k, exponential of the first protein
lb(totNumParams+1:totNumParams+1)=0; %background parameter b

ub(1:totNumParams/2)=upperBoundkFactor;
ub(totNumParams/2+1:totNumParams)=0;
ub(totNumParams+1:totNumParams+1)=1e-8;
numBoots=500; %nubmer of bootstraps estiamtes to fit on
numPreBoots=500;%Number of bootstraps used when calculating the weights for the fit 
perc=0;
tic;
[expFitCombDiff] = fitAutoCorrModelv9_bgguess(timeWind,normNumber,kguessExpMax,kguessExpMin,numGuesses,theScales,compare,lb,ub,numBoots,numPreBoots,perc,scaleRes,combUnNorm,bgInds,aguessExpMax,aguessExpMin);
%[expFitCombDiff] = fitAutoCorrModelv2(timeWind,normNumber,kguessExpMax,kguessExpMin,numGuesses,theScales,compare,lb,ub,numBoots,numPreBoots,perc,scaleRes);
toc
ind1=1;
kfacInds=1:totNumParams/2;
kInds=totNumParams/2+1:totNumParams;

bExtraInds=(totNumParams+1)*ones(1,numScales);
bootInd=8;
for i=1:numScales
    figure('position',[0,0,250,500]);
 
    hold on;
    currAutoCorr=expFitCombDiff{ind1}.autoCorrsDiff(i,:);
  
    
    currAutoCorrBoot=expFitCombDiff{ind1}.boot.autoCorrsDiffBoot(i,:,:);
    errors=std(currAutoCorrBoot,[],3);
    errors=reshape(errors,[1 timeWind]);
    errorbar(expFitCombDiff{ind1}.lagTimes(i,:),currAutoCorr,errors,'kx','LineWidth',2);
    bgAmpl=expFitCombDiff{ind1}.bgAmpl;

   hold on;
   currbestparams=expFitCombDiff{ind1}.bestparams;
   fac=10;
   timeDiff=expFitCombDiff{ind1}.lagTimes(i,2)-expFitCombDiff{ind1}.lagTimes(i,1);
   plotTimes=expFitCombDiff{ind1}.lagTimes(i,1):timeDiff/fac:expFitCombDiff{ind1}.lagTimes(i,end);
   normInds=zeros(1,normNumber);
   for j=1:normNumber
       [~,normInds(j)]=min(abs(plotTimes-expFitCombDiff{ind1}.lagTimes(i,j)));
       
   end
   autoCorrScaledBg1=combUnNorm{bgInds(1)}.autoCorrScaled;
   autoCorrScaledBg2=combUnNorm{bgInds(2)}.autoCorrScaled;
   
   
   
   bg1=[autoCorrScaledBg1(:,1)-autoCorrScaledBg1(:,end)];
   bg2=[autoCorrScaledBg2(:,1)-autoCorrScaledBg2(:,end)];
   corr1=corrFuncMultExpBgbExtra4_normIndscleaned(ones(1,numTerms),currbestparams(kInds),currbestparams(bExtraInds(i))*bg1(i),plotTimes,normInds);
   corr2=corrFuncMultExpBgbExtra4_normIndscleaned(ones(1,numTerms),currbestparams(kInds).*currbestparams(kfacInds),currbestparams(bExtraInds(i))*bg2(i),plotTimes,normInds);
   
   plot(plotTimes,corr1-corr2,'k','LineWidth',2);
    ylim([-0.3 0.3])
    xlim([0 expFitCombDiff{ind1}.lagTimes(i,end)]);
    %xlim([10^-6 10^-3]);
    xlabel('Time lag [s]');
    ylabel('Autocorrelation diff.')
    set(gca,'Xscale','log');

end


i=1; %index of timescale in theScales to do more plotting for
figure('position',[0,0,250,500]);
%figure()
%title(['Time scale: ' num2str(timeRegimes(theScales(i)))])
hold on;
currAutoCorr=expFitCombDiff{ind1}.autoCorrsDiff(i,:);
%currAutoCorr=expFitCombDiff{ind1}.boot.autoCorrsDiffBoot(i,:,bootInd);

currAutoCorrBoot=expFitCombDiff{ind1}.boot.autoCorrsDiffBoot(i,:,:);

currAutoCorrDiffCorrected=expFitCombDiff{ind1}.correctedCorr1(i,:)-expFitCombDiff{ind1}.correctedCorr2(i,:);
currAutoCorrDiffCorrectedBoot=expFitCombDiff{ind1}.boot.correctedCorr1Boot(i,:,:)-expFitCombDiff{ind1}.boot.correctedCorr2Boot(i,:,:);


errors=std(currAutoCorrBoot,[],3);
errors=reshape(errors,[1 timeWind]);
errorbar(expFitCombDiff{ind1}.lagTimes(i,:),currAutoCorr,errors,'kx','LineWidth',2);
bgAmpl=expFitCombDiff{ind1}.bgAmpl;

hold on;
currbestparams=expFitCombDiff{ind1}.bestparams;
fac=1000;
timeDiff=expFitCombDiff{ind1}.lagTimes(i,2)-expFitCombDiff{ind1}.lagTimes(i,1);
plotTimes=expFitCombDiff{ind1}.lagTimes(i,1):timeDiff/fac:expFitCombDiff{ind1}.lagTimes(i,end);
normInds=zeros(1,normNumber);
for j=1:normNumber
    [~,normInds(j)]=min(abs(plotTimes-expFitCombDiff{ind1}.lagTimes(i,j)));
    
end

pitches=[10 20 40 100];
bgs=[0 0.1 1];   
%D=nanmean(combines{compare{1}(1)}.Dxsok);

%D2=nanmean(combines{compare{1}(2)}.Dxsok);


%kFac=bestparams{1}(1);%D2/D;
kFac=kFactor;
bpFac=0.332*10^-3;%um/bp
kValues=-4*D./((pitches*bpFac/(2*pi)).^2);


signs={'-.','--','-',':'};
cols={'y','g','c','m'};
jj=1;
expFitCombDiffConstr=cell(1,numel(pitches));
for j=1:numel(pitches)
    lb(1:totNumParams/2)=lowerBoundkFactor;%kFactor*0.9999;%0
    %lb(1)=kFactor*0.9999;
    lb(totNumParams/2+1:totNumParams)=kValues(j)*1.0000001;%[-inf];%[-1e9 -1e5 -1e3];
    lb(3)=0;
    %lb(3)=0;
    
    ub(1:totNumParams/2)=upperBoundkFactor;%inf
    %ub(1)=kFactor;
    ub(totNumParams/2+1:totNumParams)=kValues(j);%[-1e5 -1e3 -1e1];
    ub(3)=1e-6;
    %ub(3)=1e-9;
    theScalesConstr=theScales(i);
    numBootsConstr=0;
    [expFitCombDiffConstr{j}] = fitAutoCorrModelv9_bgguess(timeWind,normNumber,kguessExpMax,kguessExpMin,2*numGuesses,theScalesConstr,compare,lb,ub,numBootsConstr,numPreBoots,perc,scaleRes,combUnNorm,bgInds,aguessExpMax,aguessExpMin);
    
    currbestparams=expFitCombDiffConstr{j}{ind1}.bestparams;
    
   
    bg1=[autoCorrScaledBg1(:,1)-autoCorrScaledBg1(:,end)];
    bg2=[autoCorrScaledBg2(:,1)-autoCorrScaledBg2(:,end)];
    corr1=corrFuncMultExpBgbExtra4_normInds(ones(1,numTerms),currbestparams(kInds),currbestparams(bExtraInds(i))*bg1(i),-inf,plotTimes,expFitCombDiffConstr{j}{ind1}.lagTimes(1,end),normInds);
    corr2=corrFuncMultExpBgbExtra4_normInds(ones(1,numTerms),currbestparams(kInds).*currbestparams(kfacInds),currbestparams(bExtraInds(i))*bg1(i),-inf,plotTimes,expFitCombDiffConstr{j}{ind1}.lagTimes(1,end),normInds);
    plot(plotTimes,corr1-corr2,[cols{j} signs{j}],'LineWidth',2);
end
handles=get(gca,'Children');
legend(handles(numel(pitches):-1:1),strsplit(num2str(pitches)));

%ylim([-0.15 0.15])
xlim([0 expFitCombDiff{ind1}.lagTimes(i,end)]);
%xlim([10^-6 10^-3]);
xlabel('Time lag [s]');
ylabel('Autocorrelation diff.')
set(gca,'Xscale','log');
title('Constrained fits for different pitches')
ylim([-0.3 0.3]);
xlim([0 expFitCombDiff{ind1}.lagTimes(i,end)]);

confInterval=0.67;
sorted=sort(expFitCombDiff{1}.boot.bestparamsBoot(:,2),'ascend');
minInd=round(numBoots*(1-confInterval)/2);
lowerBound=sorted(minInd);
upperBound=sorted(end-minInd);

bpFac=0.332*10^-3;%um/bp
bpMean=2*pi*sqrt(-4*D/(-10^mean(log10(-expFitCombDiff{1}.boot.bestparamsBoot(:,2)))))/bpFac;
bpLow=2*pi*sqrt(-4*D/(upperBound))/bpFac;
bpHigh=2*pi*sqrt(-4*D/(lowerBound))/bpFac;
bpReal=2*pi*sqrt(-4*D/(expFitCombDiff{1}.bestparams(:,2)))/bpFac;

sorted2=sort(expFitCombDiff{1}.boot.bestparamsBoot(:,2).*expFitCombDiff{1}.boot.bestparamsBoot(:,1),'ascend');
minInd=round(numBoots*(1-confInterval)/2);
lowerBound2=sorted2(minInd);
upperBound2=sorted2(end-minInd);
bpMean2=2*pi*sqrt(-4*D2/(-10^mean(log10(-expFitCombDiff{1}.boot.bestparamsBoot(:,2).*expFitCombDiff{1}.boot.bestparamsBoot(:,1)))))/bpFac;
bpLow2=2*pi*sqrt(-4*D2/(upperBound2))/bpFac;
bpHigh2=2*pi*sqrt(-4*D2/(lowerBound2))/bpFac;
bpReal2=2*pi*sqrt(-4*D2/(expFitCombDiff{1}.bestparams(:,2).*expFitCombDiff{1}.bestparams(:,1)))/bpFac;

theBps=2*pi*sqrt(-4*D./(expFitCombDiff{1}.boot.bestparamsBoot(:,2)))/bpFac;
theBps2=2*pi*sqrt(-4*D2./(expFitCombDiff{1}.boot.bestparamsBoot(:,2).*expFitCombDiff{1}.boot.bestparamsBoot(:,1)))/bpFac;

%Exclude the occational extreme bootstrap outlier when calculating the standard error, caused by bad fit or
%noninforamtive/close to zero autocorrelation functon in the bootstrapping
bpError=std(theBps(theBps<1e5));
bpError2=std(theBps2(theBps<1e5));


probDNAPeriod=sum(theBps<10.5)/numBoots;

disp([' The pitch of rotation coupled sliding is  : ' num2str(bpReal) ' ±' num2str(bpError) ' bp' ]);
disp([num2str(confInterval) ' confidence interval of the pitch of rotatonal coupled sliding is : [' num2str(bpHigh)  ',' num2str(bpLow) '] bp' ]);
%%
save([savePathAbs '/All_data_analysed_and_model_fitted'])