function displayHvsMMayo(folderSourceString)
close all;

if ~exist('folderSourceString','var'); folderSourceString='E:\Mayo'; end

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.8; backgroundColor = 'w'; % Panels

% Electrode Grid will now be shown according to session chosen
% Show electrodes
electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,'blank');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hParameterPanel=uipanel('Title','Parameters','FontSize',fontSizeLarge,'Position',[0.3 panelStartHeight 0.2 panelHeight]);
paramsHeight=1/3;

%fileNameString
[fileNameStringAll,fileNameStringListAll,fileNameStringListArray] = getFileNameStringList;
uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0 1-paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Session(s)','FontSize',fontSizeSmall);
hSessions=uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0.5 1-paramsHeight 0.5 paramsHeight],...
    'Style','popupmenu','String',fileNameStringAll,'fontSize',fontSizeSmall,'BackgroundColor',backgroundColor);

%neuronPopulation
populationString=[{'All'} {'Stimulated'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0 1-2*paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Population','FontSize',fontSizeSmall);
hPopulationType=uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0.5 1-2*paramsHeight 0.5 paramsHeight],...
    'Style','popupmenu','String',populationString,'fontSize',fontSizeSmall,'BackgroundColor',backgroundColor);

% Analysis Interval
analysisIntervalString = [{'Baseline'} {'StimOnset'} {'TargetOnset_250ms'} {'TargetOnset_500ms'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0 1-3*paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);
hTimePeriodType=uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0.5 1-3*paramsHeight 0.5 paramsHeight],...
    'Style','popupmenu','String',analysisIntervalString,'fontSize',fontSizeSmall,'BackgroundColor',backgroundColor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Analysis Option %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysisHeight=1/6;
tapers=[2 3];
performInterval=[50 20];
hAnalysisPanel=uipanel('Title','Analysis Options','FontSize',fontSizeLarge,'Position',[0.55 panelStartHeight 0.2 panelHeight]);
uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0 1-analysisHeight 0.5 analysisHeight],...
    'Style','text','String','MT Parameters','FontSize',fontSizeSmall);
uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.55 1-analysisHeight 0.2 analysisHeight],...
    'Style','text','String','TW','FontSize',fontSizeSmall);
uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.75 1-analysisHeight 0.2 analysisHeight],...
    'Style','text','String','K','FontSize',fontSizeSmall);
uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0 1-2*analysisHeight 0.5 analysisHeight],...
    'Style','text','String','Tapers','FontSize',fontSizeSmall);
hBandwidth=uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.55 1-2*analysisHeight 0.2 analysisHeight],...
    'Style','edit','String',num2str(tapers(1)),'FontSize',fontSizeSmall);
hTapers=uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.75 1-2*analysisHeight 0.2 analysisHeight],...
    'Style','edit','String',num2str(tapers(2)),'FontSize',fontSizeSmall);
uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0 1-3*analysisHeight 0.5 analysisHeight],...
    'Style','text','String','Performance Interval(%)','FontSize',fontSizeSmall);
hAccuracy=uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.55 1-3*analysisHeight 0.15 analysisHeight],...
    'Style','edit','String',num2str(performInterval(1)),'FontSize',fontSizeSmall);
uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.7 1-3*analysisHeight 0.1 analysisHeight],...
    'Style','text','String','+/-','FontSize',fontSizeMedium);
hTolerence=uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0.8 1-3*analysisHeight 0.15 analysisHeight],...
    'Style','edit','String',num2str(performInterval(2)),'FontSize',fontSizeSmall);
hLDA=uicontrol('Parent',hAnalysisPanel,'Units','normalized','Position',[0 1-4*analysisHeight 1 analysisHeight],...
    'Style','togglebutton','String','Linear Discriminant Analysis','FontSize',fontSizeMedium);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% plot options %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlotPanel=uipanel('Title','Plotting Options','FontSize',fontSizeLarge,'Position',[0.75 panelStartHeight 0.2 panelHeight]);

uicontrol('Parent',hPlotPanel,'Unit','Normalized','Position',[0 0.5 1 0.2],...
    'Style','pushbutton','String','Plot/Get Data(LDA)','fontSize',fontSizeMedium,'Callback',{@plotData_Callback});
uicontrol('Parent',hPlotPanel,'Unit','Normalized','Position',[0,0.75 1 0.2],...
    'Style','pushbutton','String','Cla','fontSize',fontSizeMedium,'Callback',{@cla_Callback});

hPSD=getPlotHandles(3,1,[0.05 0.05 0.2 0.7],0.05); linkaxes(hPSD);
hPSTH=getPlotHandles(3,1,[0.3 0.05 0.2 0.7],0.05); linkaxes(hPSTH);
hBarFR=getPlotHandles(3,1,[0.55 0.05 0.1 0.7],0.05); linkaxes(hBarFR);
hBarAlpha=getPlotHandles(3,1,[0.7 0.05 0.1 0.7],0.05); linkaxes(hBarAlpha);
hBarGamma=getPlotHandles(3,1,[0.85 0.05 0.1 0.7],0.05); linkaxes(hBarGamma);

    function plotData_Callback(~,~)
        
        s=get(hSessions,'val'); fileNameString = fileNameStringListArray{s} ; SessionIDString = fileNameStringListAll{s};
        p=get(hPopulationType,'val'); populationType=populationString{p};
        tp=get(hTimePeriodType,'val'); tpStr= analysisIntervalString{tp};
        accu=get(hAccuracy,'String'); accuracy= str2double(accu);
        tol=get(hTolerence,'String'); tolerence= str2double(tol);
        performInterval=[accuracy tolerence];
        LDAFlag=get(hLDA,'val');
        
        bandwidth = get(hBandwidth,'String'); TW = str2double(bandwidth);
        tapersSet = get(hTapers,'String'); K = str2double(tapersSet);
        if K~=2*TW-1
            error('Set TW and K as K = 2TW-1')
        elseif K==2*TW-1
            tapers = [TW K];
        end
        
        % Show electrodes
        if strcmp(SessionIDString{1},'all (N=24)')
            hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,'blank');
        else
            electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
            hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,SessionIDString{1});
        end
        
        [PSDDataHit,PSDDataMiss,FRHit,FRMiss,freqValsMT,electrodeArray]=analyseAllSession(fileNameString,folderSourceString,tpStr,populationType,tapers,performInterval,LDAFlag);
        if LDAFlag==0
            colorNames='bmg';
            alphaRangeHz = [8 12]; gammaRangeHz = [40 80]; lineNoiseFreqHz= 60;
            alphaPos = intersect(find(freqValsMT>=alphaRangeHz(1)),find(freqValsMT<=alphaRangeHz(2)));
            lineNoisePos = find(freqValsMT==lineNoiseFreqHz);
            gammaPos = setdiff(intersect(find(freqValsMT>=gammaRangeHz(1)),find(freqValsMT<=gammaRangeHz(2))),lineNoisePos);
            for i=1:3
                if s==length(fileNameStringListAll)||i==3
                else
                    showElectrodeLocationsMayo([],electrodeArray{i},colorNames(i),hElectrodes,1,0,SessionIDString{1}); % Show electrodes used for analysis
                end
                for cue=1:3
                    clear DataPSD
                    DataPSD=squeeze(10*(log10(PSDDataHit{i}(cue,:,:))-log10(PSDDataMiss{i}(cue,:,:))));
                    plotData(hPSD(i),freqValsMT,DataPSD,colorNames(cue));
                    line([0 100],[0 0],'Color','k','Parent',hPSD(i));
                end
                ChangeFR{i}=FRHit{i}-FRMiss{i};
                alphadB{i}=10*(log10(sum(PSDDataHit{i}(:,:,alphaPos),3))-log10(sum(PSDDataMiss{i}(:,:,alphaPos),3)));
                gammadB{i}=10*(log10(sum(PSDDataHit{i}(:,:,gammaPos),3))-log10(sum(PSDDataMiss{i}(:,:,gammaPos),3)));
                getBarPlot(hBarFR(i),ChangeFR{i},colorNames);
                getBarPlot(hBarAlpha(i),alphadB{i},colorNames);
                getBarPlot(hBarGamma(i),gammadB{i},colorNames);
            end
            
            yLims=getYLims(hPSD); axis(hPSD(1),[0 100 yLims]);
            yLims=getYLims(hBarFR); axis(hBarFR(1),[0 4 yLims]);
            yLims=getYLims(hBarAlpha); axis(hBarAlpha(1),[0 4 yLims]);
            yLims=getYLims(hBarGamma); axis(hBarGamma(1),[0 4 yLims]);
            
            %Title
            title(hPSD(1),'Power Spectra');
            title(hBarFR(1),'Firing rate');
            title(hBarAlpha(1),'Alpha');
            title(hBarGamma(1),'Gamma');
            
            %X Label
            xlabel(hPSD(3),'Frequency (Hz)');
            
            %Y Label
            arraySideString=[{'Right'} {'Left'} {'Both'}];
            for side=1:2
                ylabel(hPSD(side),[arraySideString{side} ' array' ' N=' num2str(size(FRHit{side},2))],'Color',colorNames(side));
            end
            
            ylabel(hPSD(3),{[arraySideString{3} ' array' ' N=' num2str(size(FRHit{3},2))];'Change in Power (dB)'});
            ylabel(hBarFR(3),'Change in firing rate (spikes/s)');
            ylabel(hBarAlpha(3),'Change in Power (dB)');
            ylabel(hBarGamma(3),'Change in Power (dB)');
            
            %legend
            cueString=[{'Cued'} {'Uncued'} {'Neutral'}];
            for m=1:3
                text('Parent',hPSD(1),'Units','normalized','Position',[0.7 0.35-0.1*(m-1)],'String',cueString{m},'FontSize',fontSizeMedium,'Color',colorNames(m));
            end
        end
    end
        function cla_Callback(~,~)
            claGivenPlotHandle(hPSD);
            claGivenPlotHandle(hBarFR);
            claGivenPlotHandle(hBarAlpha);
            claGivenPlotHandle(hBarGamma);
            claGivenPlotHandle(hElectrodes);
        end
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for j=1:numRows
                for k=1:numCols
                    cla(plotHandles(j,k));
                end
            end
        end
end


function [PSDDataHit,PSDDataMiss,FRHit,FRMiss,freqValsMT,electrodeArray]=analyseAllSession(fileNameString,folderSourceString,tpStr,populationType,tapers,performInterval,LDAFlag)

disp(['Working on dataset 1 of ' num2str(length(fileNameString))]);
[~,~,~,~,oriChangeIndex]=getGoodPos(fileNameString{1},folderSourceString,performInterval);
if nnz(isnan(oriChangeIndex))~=0
    disp(['Performance in  ' fileNameString{1} 'is out of the interval ' num2str(performInterval(1)) ' +/-' num2str(performInterval(2)) '%']);
else
    [PSDDataHit,PSDDataMiss,FRHit,FRMiss,freqValsMT,electrodeArray]=analyseSingleSession(fileNameString{1},folderSourceString,tpStr,populationType,tapers,performInterval,LDAFlag); 
end
if length(fileNameString)>1
    for i=2:length(fileNameString)
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameString))]);
        [~,~,~,~,oriChangeIndex]=getGoodPos(fileNameString{i},folderSourceString,performInterval);
        if nnz(isnan(oriChangeIndex))~=0
            disp(['Performance in  ' fileNameString{1} ' is out of the interval ' num2str(performInterval(1)) ' +/- ' num2str(performInterval(2)) '%']);
            continue
        else
            [PSDDataHitTMP,PSDDataMissTMP,FRHitTMP,FRMissTMP,freqValsMT,electrodeArrayTMP]=analyseSingleSession(fileNameString{i},folderSourceString,tpStr,populationType,tapers,performInterval,LDAFlag);  
        end
        if ~exist('PSDDataHit','var')
            PSDDataHit=PSDDataHitTMP;
            PSDDataMiss=PSDDataMissTMP;
            FRHit=FRHitTMP;
            FRMiss=FRMissTMP;
            electrodeArray=electrodeArrayTMP;
            continue
        end
        for side=1:3
            PSDDataHit{side}=cat(2,PSDDataHit{side},PSDDataHitTMP{side});
            PSDDataMiss{side}=cat(2,PSDDataMiss{side},PSDDataMissTMP{side});
            FRHit{side}=cat(2,FRHit{side},FRHitTMP{side});
            FRMiss{side}=cat(2,FRMiss{side},FRMissTMP{side});
        end
        for m=1:2
            electrodeArray{m}=cat(1,electrodeArray{m},electrodeArrayTMP{m});
        end
    end
end           
end

function [PSDDataHit,PSDDataMiss,FRHit,FRMiss,freqValsMT,electrodeArray]=analyseSingleSession(fileNameString,folderSourceString,tpStr,populationType,tapers,performInterval,LDAFlag)

% load Data
folderName=fullfile(folderSourceString,'Data','segmentedData',fileNameString);
attCueList=[{'0V'} {'1V'} {'0I'} {'1I'} {'N'}];

if strcmp(tpStr,'Baseline')
    timeRange=[-0.25 0];
    
    for cue=1:5
        lfpDataHit{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{cue} '_StimOnset_LFP']));
        lfpDataMiss{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{cue} '_StimOnset_LFP']));
        spikeDataHit{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{cue} '_StimOnset_Spikes']));
        spikeDataMiss{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{cue} '_StimOnset_Spikes']));
    end
    
    
elseif strcmp(tpStr,'StimOnset')
    timeRange=[0 0.25];
    
    for cue=1:5
        lfpDataHit{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{cue} '_StimOnset_LFP']));
        lfpDataMiss{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{cue} '_StimOnset_LFP']));
        spikeDataHit{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{cue} '_StimOnset_Spikes']));
        spikeDataMiss{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{cue} '_StimOnset_Spikes']));
    end
    
    
elseif strcmp(tpStr(1:11),'TargetOnset')
    if strcmp(tpStr(13:end),'500ms')
        timeRange=[-0.5 0];
    elseif strcmp(tpStr(13:end),'250ms')
        timeRange=[-0.25 0];
    end
    
    for cue=1:5
        lfpDataHit{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{cue} '_TargetOnset_LFP']));
        lfpDataMiss{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{cue} '_TargetOnset_LFP']));
        spikeDataHit{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{cue} '_TargetOnset_Spikes']));
        spikeDataMiss{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{cue} '_TargetOnset_Spikes']));
    end
    
end

timeVals=lfpDataHit{1}.timeVals;
Fs              = round(1/(timeVals(2)-timeVals(1)));
timePos= intersect(find(timeVals>=timeRange(1)),find(timeVals<timeRange(2)));
freqVals = 0:1/(diff(timeRange)):Fs-1/(diff(timeRange));
alphaRangeHz = [8 12]; gammaRangeHz = [40 80]; lineNoiseFreqHz= 60;
alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
lineNoisePos = find(freqVals==lineNoiseFreqHz);
gammaPos = setdiff(intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<=gammaRangeHz(2))),lineNoisePos);

[goodPosHit,goodPosMiss,uniqueOrientationChangeDeg,orientationChangeDeg,oriChangeIndex]= getGoodPos(fileNameString,folderSourceString,performInterval);

for cue=1:5
    oriChangeThisConditionHit=orientationChangeDeg(goodPosHit{cue});
    oriChangeThisConditionMiss=orientationChangeDeg(goodPosMiss{cue});
    selectedPosHit=find(oriChangeThisConditionHit==uniqueOrientationChangeDeg(oriChangeIndex(cue)));
    selectedPosMiss=find(oriChangeThisConditionMiss==uniqueOrientationChangeDeg(oriChangeIndex(cue)));
    lfpDataHit{cue}.segmentedLFPData=lfpDataHit{cue}.segmentedLFPData(:,selectedPosHit,timePos);
    lfpDataMiss{cue}.segmentedLFPData=lfpDataMiss{cue}.segmentedLFPData(:,selectedPosMiss,timePos);
    spikeDataHit{cue}.segmentedSpikeData=spikeDataHit{cue}.segmentedSpikeData(:,selectedPosHit);
    spikeDataMiss{cue}.segmentedSpikeData=spikeDataMiss{cue}.segmentedSpikeData(:,selectedPosMiss);
end

electrodeArray=getGoodElectrodes(fileNameString,folderSourceString,populationType);


for cue=1:3
    disp(['H' attCueList{cue} ' Stim Repeats ' num2str(size(lfpDataHit{cue}.segmentedLFPData,2))]);
    disp(['M' attCueList{cue} ' Stim Repeats ' num2str(size(lfpDataMiss{cue}.segmentedLFPData,2))]);
end

% Parameters for MT

params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 0;

folderSave = fullfile(folderSourceString,'Data','savedDataHvsM');
if ~exist(folderSave,'dir'); makeDirectory(folderSave); end
file2save=fullfile(folderSave,[fileNameString populationType tpStr '_PI' num2str(performInterval(2)) '_tapers' num2str(tapers(2))]);


for side=1:2  %1=Right  2=Left
    eList=electrodeArray{side};
    for cue=1:5  
        for k=1:length(eList)
            lfpHit=squeeze(lfpDataHit{cue}.segmentedLFPData(eList(k),:,:))';
            lfpMiss=squeeze(lfpDataMiss{cue}.segmentedLFPData(eList(k),:,:))';
            [PSDDataHit{side}{cue}(k,:,:),freqValsMT]=mtspectrumc(lfpHit,params);
            PSDDataMiss{side}{cue}(k,:,:)=mtspectrumc(lfpMiss,params);
            FRHit{side}{cue}(k,:)=getSpikeCounts(spikeDataHit{cue}.segmentedSpikeData(eList(k),:),timeRange)/diff(timeRange);
            FRMiss{side}{cue}(k,:)=getSpikeCounts(spikeDataMiss{cue}.segmentedSpikeData(eList(k),:),timeRange)/diff(timeRange);
            alphaDataHit{side}{cue}(k,:)=squeeze(sum(PSDDataHit{side}{cue}(k,alphaPos,:),2));
            alphaDataMiss{side}{cue}(k,:)=squeeze(sum(PSDDataMiss{side}{cue}(k,alphaPos,:),2));
            gammaDataHit{side}{cue}(k,:)=squeeze(sum(PSDDataHit{side}{cue}(k,gammaPos,:),2));
            gammaDataMiss{side}{cue}(k,:)=squeeze(sum(PSDDataMiss{side}{cue}(k,gammaPos,:),2));
        end
    end
end

if unique(diff(freqVals))~=unique(diff(freqValsMT))
    error('Frequency Values does not match')
end

PSDDataHit{3}=combineData(PSDDataHit);
PSDDataMiss{3}=combineData(PSDDataMiss);
FRHit{3}=combineData(FRHit);
FRMiss{3}=combineData(FRMiss);


save(file2save,'FRHit','FRMiss','alphaDataHit','alphaDataMiss','gammaDataHit','gammaDataMiss');

for side=1:2
    PSDDataHit{side}=cellfun(@(x) mean(x,3),PSDDataHit{side},'UniformOutput',false);
    PSDDataMiss{side}=cellfun(@(x) mean(x,3),PSDDataMiss{side},'UniformOutput',false);
    FRHit{side}=cellfun(@(x) mean(x,2),FRHit{side},'UniformOutput',false);
    FRMiss{side}=cellfun(@(x) mean(x,2),FRHit{side},'UniformOutput',false);
    alphaDataHit{side}=cellfun(@(x) mean(x,2),alphaDataHit{side},'UniformOutput',false);
    alphaDataMiss{side}=cellfun(@(x) mean(x,2),alphaDataMiss{side},'UniformOutput',false);
    gammaDataHit{side}=cellfun(@(x) mean(x,2),gammaDataHit{side},'UniformOutput',false);
    gammaDataMiss{side}=cellfun(@(x) mean(x,2),gammaDataMiss{side},'UniformOutput',false);
end
end


function combinedData = combineDataAcrossBothArrays(data)
if numel(size(data{1}))==2
    Data{1} = cat(2,data{1}(1,:),data{2}(2,:)); % Attend In - Valid [(R)H0V & (L)H1V]
    Data{2} = cat(2,data{1}(2,:),data{2}(1,:)); % Attend Out - Valid [(R)(H1V & (L)H0V)]
    Data{3} = cat(2,data{1}(3,:),data{2}(4,:)); % Attend In - Invalid [(R)H0I & (L)H1I]
    Data{4} = cat(2,data{1}(4,:),data{2}(3,:)); % Attend Out - Invalid [(R)H1I & (L)H0I)
    Data{5} = cat(2,data{1}(5,:),data{2}(5,:)); % Neutral    [(R)N & (L)N]
elseif numel(size(data{1}))==3
    Data{1} = cat(2,data{1}(1,:,:),data{2}(2,:,:)); % Attend In - Valid [(R)H0V & (L)H1V]
    Data{2} = cat(2,data{1}(2,:,:),data{2}(1,:,:)); % Attend Out - Valid [(R)(H1V & (L)H0V)]
    Data{3} = cat(2,data{1}(3,:,:),data{2}(4,:,:)); % Attend In - Invalid [(R)H0I & (L)H1I]
    Data{4} = cat(2,data{1}(4,:,:),data{2}(3,:,:)); % Attend Out - Invalid [(R)H1I & (L)H0I)
    Data{5} = cat(2,data{1}(5,:,:),data{2}(5,:,:)); % Neutral    [(R)N & (L)N]
end

combinedData = cat(1,Data{1},Data{2},Data{3},Data{4},Data{5});
end


function [goodPosHit,goodPosMiss,uniqueOrientationChangeDeg,orientationChangeDeg,oriChangeIndex]= getGoodPos(fileNameString,folderSourceString,performInterval)

[perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString);

goodPosHit=[goodIndexList(1:4) goodIndexList(9)];
goodPosMiss=[goodIndexList(5:8) goodIndexList(10)];


% behavior{1}=perCorrect([1 3 5],:);
% behavior{2}=perCorrect([2 4 5],:);
% behavior{3}=(behavior{1}+behavior{2})/2;
performInterval=performInterval/100;

for cue=1:5
    [performance(cue),oriChangeIndex(cue)]= min(abs(perCorrect(cue,:)-performInterval(1)));
    if performance(cue)>performInterval(2); oriChangeIndex(cue)=NaN; end
end

end
function [fileNameStringAll,fileNameStringListAll,fileNameStringListArray] = getFileNameStringList

[tmpFileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

fileNameStringAll = ''; pos=1;
clear fileNameStringListArray

for i=1:length(monkeyNameList)
    for j=1:length(tmpFileNameStringList{i})
        fileNameStringAll = [cat(2,fileNameStringAll,tmpFileNameStringList{i}{j}) '|'];
        fileNameStringListAll{pos} = tmpFileNameStringList{i}(j);
        fileNameStringListArray{pos} = tmpFileNameStringList{i}(j); %#ok<*AGROW>
        pos=pos+1;
    end
end

allNames = [];
for i=1:length(monkeyNameList)
    fileNameStringAll = [cat(2,fileNameStringAll,monkeyNameList{i}) ' (N=' num2str(length(tmpFileNameStringList{i})) ')|'];
    fileNameStringListAll{pos} = {[monkeyNameList{i} ' (N=' num2str(length(tmpFileNameStringList{i})) ')']};
    fileNameStringListArray{pos} = tmpFileNameStringList{i};
    allNames = cat(2,allNames,tmpFileNameStringList{i});
    pos=pos+1;
end

fileNameStringAll = cat(2,fileNameStringAll,['all (N=' num2str(length(allNames)) ')']);
fileNameStringListAll{pos} = {['all (N=' num2str(length(allNames)) ')']};
fileNameStringListArray{pos} = allNames;
end
function electrodeArray = getGoodElectrodes(fileNameString,folderSourceString,populationType)

suaCutoff = 3;
% Get sorting rating
sortRatings=sortRating(fileNameString,folderSourceString);
sua=intersect(find(sortRatings>0),find(sortRatings<=suaCutoff));
mua=find(sortRatings>suaCutoff);
all=union(sua,mua,'sorted');

[~,~,electrodeArrayPos]=electrodePositionOnGridMayo(1,fileNameString);

spkData{1}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H1V_StimOnset_Spikes']));
spkData{2}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H0V_StimOnset_Spikes']));

electrodeArray{1}=intersect(all,electrodeArrayPos(:,8:13));
electrodeArray{2}=intersect(all,electrodeArrayPos(:,1:6));

if strcmp(populationType,'Stimulated')
    for j=1:2 %for each array
        elecList=electrodeArray{j}(:);
        for k=1:length(elecList)
            FRBl(j,k) =mean(getSpikeCounts(spkData{j}.segmentedSpikeData(elecList(k),:),[-0.25 0]))/diff([-0.25 0]);
            FRSt(j,k) =mean(getSpikeCounts(spkData{j}.segmentedSpikeData(elecList(k),:),[0.25 0.5]))/diff([0.25 0.5]);
        end
        delFR{j}=FRSt(j,:)-FRBl(j,:);
        electrodeArray{j}=elecList(find(delFR{j}>5)); %#ok<FNDSB> % Threshold for choosing stimulated units is 5 spikes/s
    end
end
end
function sortRating=sortRating(fileNameString,folderSourceString)

[~,~,ratings]=xlsread(fullfile(folderSourceString,'Data','extractedData','SortRatings_CUonly'),fileNameString,'A1:A96');
if length(ratings)~=96
    error('electrode numbers do not match')
end

for i=1:96
    if ischar(ratings{i})
        ratings{i}=NaN;       % to eliminate electrodes rated as x, xx and blank.
    end
end
sortRating=cell2mat(ratings);
end
function getBarPlot(h,data,colorNames)
mData=mean(data,2);
semData=std(data,[],2)/sqrt(size(data,2));
for cue=1:3
    bar(h,cue,mData(cue),colorNames(cue));
    hold(h,'on');
    errorbar(h,cue,mData(cue),semData(cue),colorNames(cue))
end
set(h,'XTick',[],'XTickLabel',[]);
end
function plotData(h,x,data,colorName)
mData=mean(data,1);
semData=std(data,[],1)/sqrt(size(data,1));
xsLong = [x fliplr(x)];
ysLong = [mData+semData fliplr(mData-semData)];

plot(h,x,mData,'color',colorName);
patch(xsLong,ysLong,colorName,'FaceAlpha',0.3,'Parent',h);
hold(h,'on')
end
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;
for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end