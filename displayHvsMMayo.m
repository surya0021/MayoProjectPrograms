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
paramsHeight=1/4;

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


%orientationChange
oriChangeString=[{'Single'} {'Selected'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0 1-3*paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Orientation','FontSize',fontSizeSmall);
hOrientation=uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0.5 1-3*paramsHeight 0.5 paramsHeight],...
    'Style','popupmenu','String',oriChangeString,'fontSize',fontSizeSmall,'BackgroundColor',backgroundColor,'Callback',@enableDisable_Callback);


% Analysis Interval
analysisIntervalString = [{'Baseline'} {'StimOnset'} {'TargetOnset_250ms'} {'TargetOnset_500ms'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0 1-4*paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);
hTimePeriodType=uicontrol('Parent',hParameterPanel,'Unit','Normalized','Position',[0.5 1-4*paramsHeight 0.5 paramsHeight],...
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
hNeutralTrials=uicontrol('Parent',hAnalysisPanel,'Unit','Normalized','Position',[0 1-4*analysisHeight 1 analysisHeight],...
    'Style','togglebutton','String','Divide Neutral Trials depending on Target Location','fontSize',fontSizeSmall,'Callback',@plotHandle_Callback);

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% plot options %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlotPanel=uipanel('Title','Plotting Options','FontSize',fontSizeLarge,'Position',[0.75 panelStartHeight 0.2 panelHeight]);

uicontrol('Parent',hPlotPanel,'Unit','Normalized','Position',[0.05 0.5 0.2 0.2],...
    'Style','pushbutton','String','Plot','fontSize',fontSizeMedium,'Callback',{@plotData_Callback});
uicontrol('Parent',hPlotPanel,'Unit','Normalized','Position',[0.05,0.75 0.2 0.2],...
    'Style','pushbutton','String','Cla','fontSize',fontSizeMedium,'Callback',{@cla_Callback});
% Figures
figure(1) % Plots of differences between hits and miss trials
hDelPSD=getPlotHandles(3,1,[0.05 0.05 0.2 0.7],0.05); linkaxes(hDelPSD);
hDelPSTH=getPlotHandles(3,1,[0.3 0.05 0.2 0.7],0.05); linkaxes(hDelPSTH);
hBarDelFR=getPlotHandles(3,1,[0.55 0.05 0.1 0.7],0.05); linkaxes(hBarDelFR);
hBarDelAlpha=getPlotHandles(3,1,[0.7 0.05 0.1 0.7],0.05); linkaxes(hBarDelAlpha);
hBarDelGamma=getPlotHandles(3,1,[0.85 0.05 0.1 0.7],0.05); linkaxes(hBarDelGamma);

fig2=figure(2); %figure 2 for Plots of Absolute values of Hits and Miss Trials
set(fig2,'units','normalized','outerposition',[0 0 1 1]); 
hPSD=getPlotHandles(3,3,[0.05 0.05 0.4 0.9],0.01,0.1);
hPSTH=getPlotHandles(3,3,[0.5 0.05 0.48 0.9],0.01,0.1);

    

    function plotData_Callback(~,~)
        
        s=get(hSessions,'val'); fileNameString = fileNameStringListArray{s} ; SessionIDString = fileNameStringListAll{s};
        p=get(hPopulationType,'val'); populationType=populationString{p};
        tp=get(hTimePeriodType,'val'); tpStr= analysisIntervalString{tp};
        o=get(hOrientation,'val');  oStr=oriChangeString{o};
        accu=get(hAccuracy,'String'); accuracy= str2double(accu);
        tol=get(hTolerence,'String'); tolerence= str2double(tol);
        performInterval=[accuracy tolerence];
        neutralFlag=get(hNeutralTrials,'val');
        
        bandwidth = get(hBandwidth,'String'); TW = str2double(bandwidth);
        tapersSet = get(hTapers,'String'); K = str2double(tapersSet);
        if K~=2*TW-1
            error('Set TW and K as K = 2TW-1')
        elseif K==2*TW-1
            tapers = [TW K];
        end
        
        % Show electrodes
        if strcmp(SessionIDString{1}(1:3),'all')
            hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,'blank');
        else
            electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
            hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,SessionIDString{1});
        end
        
        [PSDDataHit,PSDDataMiss,psthHit,psthMiss,FRHit,FRMiss,freqVals,xsFR,electrodeArray]=analyseAllSession(fileNameString,folderSourceString,tpStr,oStr,populationType,tapers,performInterval,neutralFlag);
        if neutralFlag
            colorNames='bmgc';
        else
            colorNames='bmg';
        end
        alphaRangeHz = [8 12]; gammaRangeHz = [40 80]; lineNoiseFreqHz= 60;
        alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
        lineNoisePos = find(freqVals==lineNoiseFreqHz);
        gammaPos = setdiff(intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<=gammaRangeHz(2))),lineNoisePos);
        for i=1:3
            if s==length(fileNameStringListAll)||i==3
            else
                showElectrodeLocationsMayo([],electrodeArray{i},colorNames(i),hElectrodes,1,0,SessionIDString{1}); % Show electrodes used for analysis
            end
            for cue=1:size(PSDDataHit{i},1) 
                clear DataPSD DataPSTH
                DataPSD=squeeze(10*(log10(PSDDataHit{i}(cue,:,:))-log10(PSDDataMiss{i}(cue,:,:))));
                DataPSTH=squeeze(psthHit{i}(cue,:,:)-psthMiss{i}(cue,:,:));
                plotData(hDelPSD(i),freqVals,DataPSD,colorNames(cue));
                plotData(hDelPSTH(i),xsFR,DataPSTH,colorNames(cue));
                line([0 100],[0 0],'Color','k','Parent',hDelPSD(i));
                line([xsFR(1) xsFR(end)],[0 0],'Color','k','Parent',hDelPSTH(i));
            end
            ChangeFR{i}=FRHit{i}-FRMiss{i};
            alphadB{i}=10*(log10(sum(PSDDataHit{i}(:,:,alphaPos),3))-log10(sum(PSDDataMiss{i}(:,:,alphaPos),3)));
            gammadB{i}=10*(log10(sum(PSDDataHit{i}(:,:,gammaPos),3))-log10(sum(PSDDataMiss{i}(:,:,gammaPos),3)));
            getBarPlot(hBarDelFR(i),ChangeFR{i},colorNames);
            getBarPlot(hBarDelAlpha(i),alphadB{i},colorNames);
            getBarPlot(hBarDelGamma(i),gammadB{i},colorNames);
        end
        
        for cue=1:size(PSDDataHit{1},1)
            for side=1:3
                plotHvsM(hPSD(3*(cue-1)+side),freqVals,log10(squeeze(PSDDataHit{side}(cue,:,:))),log10(squeeze(PSDDataMiss{side}(cue,:,:))));
                plotHvsM(hPSTH(3*(cue-1)+side),xsFR,squeeze(psthHit{side}(cue,:,:)),squeeze(psthMiss{side}(cue,:,:)));
            end
        end
        yLims=getYLims(hDelPSD); axis(hDelPSD(1),[0 100 yLims]);
        yLims=getYLims(hBarDelFR); axis(hBarDelFR(1),[0 5 yLims]);
        yLims=getYLims(hBarDelAlpha); axis(hBarDelAlpha(1),[0 5 yLims]);
        yLims=getYLims(hBarDelGamma); axis(hBarDelGamma(1),[0 5 yLims]);
        
        %Title
        title(hDelPSD(1),'Power Spectra');
        title(hBarDelFR(1),'Firing rate');
        title(hBarDelAlpha(1),'Alpha');
        title(hBarDelGamma(1),'Gamma');
        title(hPSD(1),'Valid');
        title(hPSTH(1),'Valid');
        title(hPSD(4),'Invalid');
        title(hPSTH(4),'Invalid');
        if size(hPSD,2)==3
            title(hPSD(7),'Neutral');
            title(hPSTH(7),'Neutral');
        elseif size(hPSD,2)==4
            title(hPSD(7),'Neutral(Contra-Change)');
            title(hPSTH(7),'Neutral(Contra-Change)');
            title(hPSD(10),'Neutral(Ipsi-Change)');
            title(hPSTH(10),'Neutral(Ipsi-Change)');
        end
        
        %X Label
        xlabel(hDelPSD(3),'Frequency (Hz)');
        xlabel(hDelPSTH(3),'Time(s)');
        xlabel(hPSD(3),'Frequency (Hz)');
        xlabel(hPSTH(3),'Time(s)');
        
        
        %Y Label
        arraySideString=[{'Right'} {'Left'} {'Both'}];
        for side=1:2
            ylabel(hDelPSD(side),[arraySideString{side} ' array' ' N=' num2str(size(FRHit{side},2))],'Color',colorNames(side));
            ylabel(hPSD(side),[arraySideString{side} ' array' ' N=' num2str(size(FRHit{side},2))],'Color',colorNames(side));
        end
        
        ylabel(hDelPSD(3),{[arraySideString{3} ' array' ' N=' num2str(size(FRHit{3},2))];'Change in Power (dB)'});
        ylabel(hPSD(3),{[arraySideString{3} ' array' ' N=' num2str(size(FRHit{3},2))];'log_{10}(Power)'});
        ylabel(hBarDelFR(3),'Change in firing rate (spikes/s)');
        ylabel(hBarDelAlpha(3),'Change in Power (dB)');
        ylabel(hBarDelGamma(3),'Change in Power (dB)');
        ylabel(hPSTH(3),'Firing Rate (spikes/s)');
        
        %legend
        if neutralFlag
            cueString=[{'Cued'} {'Uncued'} {'Neutral(Contra)'} {'Neutral(Ipsi)'}];
        else
            cueString=[{'Cued'} {'Uncued'} {'Neutral'}];
        end
        for m=1:length(cueString)
            text('Parent',hDelPSD(1),'Units','normalized','Position',[0.7 0.45-0.1*(m-1)],'String',cueString{m},'FontSize',fontSizeSmall,'Color',colorNames(m));
        end
        text('Parent',hPSD(1),'Units','normalized','Position',[0.7 0.85],'String','Hit','FontSize',fontSizeMedium,'Color','r');
        text('Parent',hPSD(1),'Units','normalized','Position',[0.7 0.75],'String','Miss','FontSize',fontSizeMedium,'Color','b');
        text('Parent',hPSTH(1),'Units','normalized','Position',[0.7 0.85],'String','Hit','FontSize',fontSizeMedium,'Color','r');
        text('Parent',hPSTH(1),'Units','normalized','Position',[0.7 0.75],'String','Miss','FontSize',fontSizeMedium,'Color','b');
    end
    function cla_Callback(~,~)
        claGivenPlotHandle(hDelPSD);
        claGivenPlotHandle(hDelPSTH);
        claGivenPlotHandle(hBarDelFR);
        claGivenPlotHandle(hBarDelAlpha);
        claGivenPlotHandle(hBarDelGamma);
        claGivenPlotHandle(hPSD);
        claGivenPlotHandle(hPSTH);
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
    function enableDisable_Callback(~,~)
        if strcmp({'Selected'},oriChangeString(hOrientation.Value));  hAccuracy.Enable='off';   hTolerence.Enable='off';
        else hAccuracy.Enable='on';  hTolerence.Enable='on';
        end
    end

    function plotHandle_Callback(~,~)
        if get(hNeutralTrials,'val')==1
            clf(fig2)
            set(0,'CurrentFigure',fig2)
            hPSD=getPlotHandles(3,4,[0.05 0.05 0.4 0.9],0.01,0.1);
            hPSTH=getPlotHandles(3,4,[0.5 0.05 0.48 0.9],0.01,0.1);
        else
           clf(fig2)
           set(0,'CurrentFigure',fig2) 
           hPSD=getPlotHandles(3,3,[0.05 0.05 0.4 0.9],0.01,0.1);
           hPSTH=getPlotHandles(3,3,[0.5 0.05 0.48 0.9],0.01,0.1);
        end
    end
end


function [PSDDataHit,PSDDataMiss,psthHit,psthMiss,FRHit,FRMiss,freqVals,xsFR,electrodeArray]=analyseAllSession(fileNameString,folderSourceString,tpStr,oStr,populationType,tapers,performInterval,neutralFlag)

disp(['Working on dataset 1 of ' num2str(length(fileNameString)) ': ' fileNameString{1}]);
if strcmp(oStr,'Single')
    [~,~,~,~,~,oriChangeIndex]=getGoodPos(fileNameString{1},folderSourceString,performInterval,neutralFlag);
else
    oriChangeIndex=[];
end
if nnz(isnan(oriChangeIndex))~=0
    disp(['Performance in  ' fileNameString{1} 'is out of the interval ' num2str(performInterval(1)) ' +/-' num2str(performInterval(2)) '%']);
else
    [PSDDataHit,PSDDataMiss,psthHit,psthMiss,FRHit,FRMiss,freqVals,xsFR,electrodeArray]=analyseSingleSession(fileNameString{1},folderSourceString,tpStr,oStr,populationType,tapers,performInterval,neutralFlag);
end

if length(fileNameString)>1
    for i=2:length(fileNameString)
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameString)) ': ' fileNameString{i}]);
        if strcmp(oStr,'Single')
            [~,~,~,~,~,oriChangeIndex]=getGoodPos(fileNameString{i},folderSourceString,performInterval,neutralFlag);
        else
            oriChangeIndex=[];
        end
        if nnz(isnan(oriChangeIndex))~=0
            disp(['Performance in  ' fileNameString{i} ' is out of the interval ' num2str(performInterval(1)) ' +/- ' num2str(performInterval(2)) '%']);
            continue
        else
            [PSDDataHitTMP,PSDDataMissTMP,psthHitTMP,psthMissTMP,FRHitTMP,FRMissTMP,freqVals,xsFR,electrodeArrayTMP]=analyseSingleSession(fileNameString{i},folderSourceString,tpStr,oStr,populationType,tapers,performInterval,neutralFlag);  
        end
        if ~exist('PSDDataHit','var')
            PSDDataHit=PSDDataHitTMP;
            PSDDataMiss=PSDDataMissTMP;
            psthHit=psthHitTMP;
            psthMiss=psthMissTMP;
            FRHit=FRHitTMP;
            FRMiss=FRMissTMP;
            electrodeArray=electrodeArrayTMP;
            continue
        end
        for side=1:3
            PSDDataHit{side}=cat(2,PSDDataHit{side},PSDDataHitTMP{side});
            PSDDataMiss{side}=cat(2,PSDDataMiss{side},PSDDataMissTMP{side});
            psthHit{side}=cat(2,psthHit{side},psthHitTMP{side});
            psthMiss{side}=cat(2,psthMiss{side},psthMissTMP{side});
            FRHit{side}=cat(2,FRHit{side},FRHitTMP{side});
            FRMiss{side}=cat(2,FRMiss{side},FRMissTMP{side});
        end
        for m=1:2
            electrodeArray{m}=cat(1,electrodeArray{m},electrodeArrayTMP{m});
        end
    end
end           
end

function [PSDDataHit,PSDDataMiss,psthHit,psthMiss,FRHit,FRMiss,freqVals,xsFR,electrodeArray]=analyseSingleSession(fileNameString,folderSourceString,tpStr,oStr,populationType,tapers,performInterval,neutralFlag)

% load Data
folderName=fullfile(folderSourceString,'Data','segmentedData',fileNameString);
if neutralFlag
    attCueList{1}=[{'0V'} {'0I'} {'0N'} {'1N'}];
    attCueList{2}=[{'1V'} {'1I'} {'1N'} {'0N'}];
else
    attCueList{1}=[{'0V'} {'0I'} {'N'}];
    attCueList{2}=[{'1V'} {'1I'} {'N'}];
end

if strcmp(tpStr,'Baseline')
    timeRange=[-0.25 0];
    for side=1:2
        for cue=1:length(attCueList{side})
            lfpDataHit{side}{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{side}{cue} '_StimOnset_LFP']));
            lfpDataMiss{side}{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{side}{cue} '_StimOnset_LFP']));
            spikeDataHit{side}{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{side}{cue} '_StimOnset_Spikes']));
            spikeDataMiss{side}{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{side}{cue} '_StimOnset_Spikes']));
        end
    end
    
elseif strcmp(tpStr,'StimOnset')
    timeRange=[0 0.25];
    for side=1:2
        for cue=1:length(attCueList{side})
            lfpDataHit{side}{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{side}{cue} '_StimOnset_LFP']));
            lfpDataMiss{side}{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{side}{cue} '_StimOnset_LFP']));
            spikeDataHit{side}{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{side}{cue} '_StimOnset_Spikes']));
            spikeDataMiss{side}{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{side}{cue} '_StimOnset_Spikes']));
        end
    end
    
elseif strcmp(tpStr(1:11),'TargetOnset')
    if strcmp(tpStr(13:end),'500ms')
        timeRange=[-0.5 0];
    elseif strcmp(tpStr(13:end),'250ms')
        timeRange=[-0.25 0];
    end
    for side=1:2
        for cue=1:length(attCueList{side})
            lfpDataHit{side}{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{side}{cue} '_TargetOnset_LFP']));
            lfpDataMiss{side}{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{side}{cue} '_TargetOnset_LFP']));
            spikeDataHit{side}{cue}= load(fullfile(folderName,[fileNameString 'H' attCueList{side}{cue} '_TargetOnset_Spikes']));
            spikeDataMiss{side}{cue}= load(fullfile(folderName,[fileNameString 'M' attCueList{side}{cue} '_TargetOnset_Spikes']));
        end
    end
end
binWidth = 10; % for PSTH
timeVals = lfpDataHit{1}{1}.timeVals;
timePos = intersect(find(timeVals>=timeRange(1)),find(timeVals<timeRange(2)));
[goodPosHit,goodPosMiss,uniqueOrientationChangeDeg,orientationChangeDeg,~,oriChangeIndex]= getGoodPos(fileNameString,folderSourceString,performInterval,neutralFlag);
 
for side=1:2
    for cue=1:length(attCueList{side})
        oriChangeThisConditionHit=orientationChangeDeg(goodPosHit{side}{cue});
        oriChangeThisConditionMiss=orientationChangeDeg(goodPosMiss{side}{cue});
        if strcmp(oStr,'Single')
            selectedPosHit=find(oriChangeThisConditionHit==uniqueOrientationChangeDeg(oriChangeIndex(side,cue)));
            selectedPosMiss=find(oriChangeThisConditionMiss==uniqueOrientationChangeDeg(oriChangeIndex(side,cue)));
        elseif strcmp(oStr,'Selected')
            selectedPosHit=sort([find(oriChangeThisConditionHit==uniqueOrientationChangeDeg(2)) find(oriChangeThisConditionHit==uniqueOrientationChangeDeg(3))]);
            selectedPosMiss=sort([find(oriChangeThisConditionMiss==uniqueOrientationChangeDeg(2)) find(oriChangeThisConditionMiss==uniqueOrientationChangeDeg(3))]);
        end
        
        lfpDataHit{side}{cue}.segmentedLFPData=lfpDataHit{side}{cue}.segmentedLFPData(:,selectedPosHit,timePos);
        lfpDataMiss{side}{cue}.segmentedLFPData=lfpDataMiss{side}{cue}.segmentedLFPData(:,selectedPosMiss,timePos);
        spikeDataHit{side}{cue}.segmentedSpikeData=spikeDataHit{side}{cue}.segmentedSpikeData(:,selectedPosHit);
        spikeDataMiss{side}{cue}.segmentedSpikeData=spikeDataMiss{side}{cue}.segmentedSpikeData(:,selectedPosMiss);
        
    end
end
electrodeArray=getGoodElectrodes(fileNameString,folderSourceString,populationType);

for side=1:2
    for cue=1:length(attCueList{side})
        disp(['H' attCueList{side}{cue} ' Stim Repeats ' num2str(size(lfpDataHit{side}{cue}.segmentedLFPData,2))]);
        disp(['M' attCueList{side}{cue} ' Stim Repeats ' num2str(size(lfpDataMiss{side}{cue}.segmentedLFPData,2))]);
    end
end
% Parameters for MT

Fs              = round(1/(timeVals(2)-timeVals(1)));
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 1;
for side=1:2    %1:Right Array  2:Left Array
    eList=electrodeArray{side};
    for cue=1:length(attCueList{side})    % For side 1, cue=[0V,0I,0N,1N] and for side 2, cue=[1V,1I,1N,0N] 
        for k=1:length(eList)
            lfpHit=squeeze(lfpDataHit{side}{cue}.segmentedLFPData(eList(k),:,:))';
            lfpMiss=squeeze(lfpDataMiss{side}{cue}.segmentedLFPData(eList(k),:,:))';
            [PSDDataHit{side}(cue,k,:),freqVals]=mtspectrumc(lfpHit,params);
            PSDDataMiss{side}(cue,k,:)=mtspectrumc(lfpMiss,params);
            [psthHit{side}(cue,k,:),xsFRHit]=getPSTH(spikeDataHit{side}{cue}.segmentedSpikeData(eList(k),:),binWidth,[timeVals(1) timeVals(end)]);
            [psthMiss{side}(cue,k,:),xsFRMiss]=getPSTH(spikeDataMiss{side}{cue}.segmentedSpikeData(eList(k),:),binWidth,[timeVals(1) timeVals(end)]);
            FRHit{side}(cue,k)=mean(getSpikeCounts(spikeDataHit{side}{cue}.segmentedSpikeData(eList(k),:),timeRange)/diff(timeRange),2);
            FRMiss{side}(cue,k)=mean(getSpikeCounts(spikeDataMiss{side}{cue}.segmentedSpikeData(eList(k),:),timeRange)/diff(timeRange),2);
        end
    end
end
if ~isequal(xsFRHit,xsFRMiss)
    error('Time points of PSTH do not match')
elseif isequal(xsFRHit,xsFRMiss)
    xsFR=xsFRHit;
end

PSDDataHit{3}=combineData(PSDDataHit,neutralFlag);
PSDDataMiss{3}=combineData(PSDDataMiss,neutralFlag);
psthHit{3}=combineData(psthHit,neutralFlag);
psthMiss{3}=combineData(psthMiss,neutralFlag);
FRHit{3}=combineData(FRHit,neutralFlag);
FRMiss{3}=combineData(FRMiss,neutralFlag);
end


function combinedData=combineData(data,neutralFlag)
% if neutralFlag
%     if numel(size(data{1}))==2
%         Data{1}=cat(2,data{1}(1:2,:),data{2}(1:2,:));   % concatenates respective Attend Valid and Attend Invalid condition of both arrays.
%         Data{2}=cat(2,data{1}(3,:),data{2}(4,:));   % concatnates contralateral changes of neutral trials i.e. left change for Right array and right change for Left Array.
%         Data{3}=cat(2,data{1}(4,:),data{2}(3,:));   % concatenates ipsilateral changes of neutral trials i.e right change for right array and left change for right array.
%     elseif numel(size(data{1}))==3
%         Data{1}=cat(2,data{1}(1:2,:,:),data{2}(1:2,:,:));   % concatenates respective Attend Valid and Attend Invalid condition of both arrays.
%         Data{2}=cat(2,data{1}(3,:,:),data{2}(4,:,:));   % concatnates contralateral changes of neutral trials i.e. left change for Right array and right change for Left Array.
%         Data{3}=cat(2,data{1}(4,:,:),data{2}(3,:,:));   % concatenates ipsilateral changes of neutral trials i.e right change for right array and left change for right array.
%         
%     end
%     combinedData= cat(1,Data{1},Data{2},Data{3});
% else
    combinedData=cat(2,data{1},data{2});
% end
end


function [goodPosHit,goodPosMiss,uniqueOrientationChangeDeg,orientationChangeDeg,behavior,oriChangeIndex]= getGoodPos(fileNameString,folderSourceString,performInterval,neutralFlag)

[perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString,neutralFlag);
if neutralFlag
    goodPosHit{1}=[goodIndexList(1:2:3) goodIndexList(9) goodIndexList(10)];    % Positioned in the order of attCueList{1}
    goodPosHit{2}=[goodIndexList(2:2:4) goodIndexList(10) goodIndexList(9)];    % Positioned in the order of attCueList{2}
    goodPosMiss{1}=[goodIndexList(5:2:7) goodIndexList(11) goodIndexList(12)];
    goodPosMiss{2}=[goodIndexList(6:2:8) goodIndexList(12) goodIndexList(11)];
else
    goodPosHit{1}=[goodIndexList(1:2:3) goodIndexList(9)];
    goodPosHit{2}=[goodIndexList(2:2:4) goodIndexList(9)];
    goodPosMiss{1}=[goodIndexList(5:2:7) goodIndexList(10)];
    goodPosMiss{2}=[goodIndexList(6:2:8) goodIndexList(10)];
end

if neutralFlag
    behavior{1}=perCorrect([1 3 5 6],:);  % 1:0V  3:0I  5:0N  6:1N  
    behavior{2}=perCorrect([2 4 6 5],:);  % 2:1V  4:1I  6:1N  5:0N
    behavior{3}=(behavior{1}+behavior{2})/2;
else
    behavior{1}=perCorrect([1 3 5],:);
    behavior{2}=perCorrect([2 4 5],:);
    behavior{3}=(behavior{1}+behavior{2})/2;
end
performInterval=performInterval/100;
for side=1:3
    for cue=1:length(goodPosHit{1})
        [performance(side,cue),oriChangeIndex(side,cue)]= min(abs(behavior{side}(cue,:)-performInterval(1)));
        if performance(side,cue)>performInterval(2); oriChangeIndex(side,cue)=NaN; end
    end
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

electrodeArray{1}=intersect(all,electrodeArrayPos(:,8:13)); % Right Array
electrodeArray{2}=intersect(all,electrodeArrayPos(:,1:6));  % left Array

if strcmp(populationType,'Stimulated')
    spkData{1}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H1V_StimOnset_Spikes']));
    spkData{2}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H0V_StimOnset_Spikes']));
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
for cue=1:size(data,1)
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

plot(h,x,mData,'color',colorName,'LineWidth',1);
patch(xsLong,ysLong,colorName,'EdgeColor','none','FaceAlpha',0.3,'Parent',h);
hold(h,'on')
end
function plotHvsM(h,x,dataHit,dataMiss)
mDataHit=mean(dataHit,1);   mDataMiss=mean(dataMiss,1);
semDataHit=std(dataHit,[],1)/sqrt(size(dataHit,1));     semDataMiss=std(dataMiss,[],1)/sqrt(size(dataMiss,1));
xsLong = [x fliplr(x)];
ysLongHit = [mDataHit+semDataHit fliplr(mDataHit-semDataHit)];      ysLongMiss = [mDataMiss+semDataMiss fliplr(mDataMiss-semDataMiss)];   

plot(h,x,mDataHit,'Color','r','LineWidth',1)
patch(xsLong,ysLongHit,'r','EdgeColor','none','FaceAlpha',0.3,'Parent',h);
hold (h,'on')
plot(h,x,mDataMiss,'Color','b','LineWidth',1)
patch(xsLong,ysLongMiss,'b','EdgeColor','none','FaceAlpha',0.3,'Parent',h);
xlim(h,[min(x) max(x)])
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