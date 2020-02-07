% Displays data from a single electrode

function displaySingleChannelGRFMayo(fileNameString,activeSide,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end

% Folders and file names
folderName = fullfile(folderSourceString,'Data','segmentedData',fileNameString);

% Get Combinations
p = load(fullfile(folderName,['parameterCombinations' num2str(activeSide) '.mat']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.3; panelStartHeight = 0.65;
staticPanelWidth = 0.25; staticStartPos = 0.025;
dynamicPanelWidth = 0.25; dynamicStartPos = 0.275;
timingPanelWidth = 0.25; timingStartPos = 0.525;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.775;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%staticTitle = [subjectName '_' expDate '_' protocolName];
% if 0 % don't plot the static panel
%     hStaticPanel = uipanel('Title','Information','fontSize', fontSizeLarge, ...
%         'Unit','Normalized','Position',[staticStartPos panelStartHeight staticPanelWidth panelHeight]);
% 
%     staticText = [{ '   '};
%         {['Monkey Name: ' subjectName]}; ...
%         {['Date: ' expDate]}; ...
%         {['Protocol Name: ' protocolName]}; ...
%         {'   '}
%         {['Orientation  (Deg): ' num2str(stimResults.orientation)]}; ...
%         {['Spatial Freq (CPD): ' num2str(stimResults.spatialFrequency)]}; ...
%         {['Eccentricity (Deg): ' num2str(stimResults.eccentricity)]}; ...
%         {['Polar angle  (Deg): ' num2str(stimResults.polarAngle)]}; ...
%         {['Sigma        (Deg): ' num2str(stimResults.sigma)]}; ...
%         {['Radius       (Deg): ' num2str(stimResults.radius)]}; ...
%         ];
% 
%     tStaticText = uicontrol('Parent',hStaticPanel,'Unit','Normalized', ...
%         'Position',[0 0 1 1], 'Style','text','String',staticText,'FontSize',fontSizeSmall);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.1; dynamicGap=0; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Analog channel
[analogChannelStringList,analogChannelStringArray] = getAnalogStringFromValues(1:96,[]);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Analog Channel','FontSize',fontSizeSmall);
hAnalogChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analogChannelStringList,'FontSize',fontSizeSmall);

% Neural channel
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Neural Channel','FontSize',fontSizeSmall);

neuralChannelsStored=[];
if ~isempty(neuralChannelsStored)
    neuralChannelString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs);
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',neuralChannelString,'FontSize',fontSizeSmall);
else
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position', [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Not found','FontSize',fontSizeSmall);
end
% Sigma
sigmaString = getStringFromValues(p.sValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeSmall);
hSigma = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeSmall);

% Spatial Frequency
spatialFreqString = getStringFromValues(p.fValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeSmall);
hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeSmall);

% Orientation
orientationString = getStringFromValues(p.oValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeSmall);
hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeSmall);

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|Raster|FFT|delta FFT|STA';
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);

% For orientation and SF
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-7*(dynamicHeight+dynamicGap) 1 dynamicHeight],...
    'Style','text','String','For sigma,ori,SF,C & TF plots','FontSize',fontSizeSmall);
% Azimuth
azimuthString = getStringFromValues(p.aValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-8*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Azimuth (Deg)','FontSize',fontSizeSmall);
hAzimuth = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-8*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',azimuthString,'FontSize',fontSizeSmall);

% Elevation
elevationString = getStringFromValues(p.eValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-9*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Elevation (Deg)','FontSize',fontSizeSmall);
hElevation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-9*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',elevationString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [0 0.4];
fftRange = [0 250];
baseline = [0 0.2];
stimPeriod = [0.05 0.25];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

% STA length
staLen = [-0.05 0.05]; 
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','STA len (s)','FontSize',fontSizeSmall);
hSTAMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(1)),'FontSize',fontSizeSmall);
hSTAMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(2)),'FontSize',fontSizeSmall);
hRemoveMeanSTA = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0 1-10*timingHeight 1 timingHeight], ...
    'Style','togglebutton','String','remove mean STA','FontSize',fontSizeMedium);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeSmall);

hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plots and message handles

% Get electrode array information
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0);

% Make plot for RFMap, centerRFMap and main Map
if length(p.aValsUnique)>=5
    mapRatio = 2/3; % this sets the relative ratio of the mapping plots versus orientation plots
else
    mapRatio = 1/2;
end

startXPos = staticStartPos; endXPos = 0.95; startYPos = 0.05; mainRFHeight = 0.55; centerGap = 0.05;
mainRFWidth = mapRatio*(endXPos-startXPos-centerGap);
otherPlotsWidth = (1-mapRatio)*(endXPos-startXPos-centerGap);

% RF and centerRF
RFMapPos = [endXPos-mainRFWidth startYPos+(3/4)*mainRFHeight mainRFWidth/2 (1/4)*mainRFHeight];
hRFMapPlot = subplot('Position',RFMapPos,'XTickLabel',[],'YTickLabel',[],'box','on');
centerRFMapPos = [endXPos-mainRFWidth+mainRFWidth/2 startYPos+(3/4)*mainRFHeight mainRFWidth/2 (1/4)*mainRFHeight];
hcenterRFMapPlot = subplot('Position',centerRFMapPos,'XTickLabel',[],'YTickLabel',[],'box','on');

% Main plot handles 
numRows = length(p.eValsUnique); numCols = length(p.aValsUnique);
gridPos=[endXPos-mainRFWidth startYPos mainRFWidth 0.65*mainRFHeight]; gap = 0.002;
plotHandles = getPlotHandles(numRows,numCols,gridPos,gap);

uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String',fileNameString,'FontSize',fontSizeLarge);

% Other functions

% Remaining Grid size
remainingWidth = otherPlotsWidth;
remainingHeight= mainRFHeight;

otherGapSize = 0.04;
otherHeight = (remainingHeight-2*otherGapSize)/3;

spatialFreqGrid  = [startXPos startYPos+ 0*(otherHeight+otherGapSize) remainingWidth otherHeight];
orientationGrid  = [startXPos startYPos+ 1*(otherHeight+otherGapSize) remainingWidth otherHeight];
sigmaGrid        = [startXPos startYPos+ 2*(otherHeight+otherGapSize) remainingWidth otherHeight];

% Plot handles
hOrientationPlot  = getPlotHandles(1,length(p.oValsUnique),orientationGrid,0.002);
hSpatialFreqPlot  = getPlotHandles(1,length(p.fValsUnique),spatialFreqGrid,0.002);
hSigmaPlot        = getPlotHandles(1,length(p.sValsUnique),sigmaGrid,0.002);

% Load analog Data
lfpData = load(fullfile(folderName,'lfpData.mat'));

goodStimNums = load(fullfile(folderName,['goodStimNums' num2str(activeSide) '.mat']));
lfpData.segmentedLFPData = lfpData.segmentedLFPData(:,goodStimNums.goodStimNums,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        a=get(hAzimuth,'val');
        e=get(hElevation,'val');
        s=get(hSigma,'val');
        f=get(hSpatialFreq,'val');
        o=get(hOrientation,'val');

        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];
        staRange = [str2double(get(hSTAMin,'String')) str2double(get(hSTAMax,'String'))];
        holdOnState = get(hHoldOn,'val');
        removeMeanSTA = get(hRemoveMeanSTA,'val');

        if analysisType==6 % Spike triggered average - do Later
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            spikeChannelPos = get(hNeuralChannel,'val');
            spikeChannelNumber = neuralChannelsStored(spikeChannelPos);
            unitID = SourceUnitIDs(spikeChannelPos);
            
            plotColors{1} = 'g';
            plotColors{2} = 'k';
            plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,blRange,stRange,folderName,staRange,removeMeanSTA);
            
            % Write code for this
            %plotSTA1Parameter1Channel(hOrientationPlot,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
            %    a,e,s,f,[],timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName);
            %plotSTA1Parameter1Channel(hSpatialFreqPlot,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
            %    a,e,s,[],o,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName);
            
            if analogChannelPos<=length(analogChannelsStored)
                analogChannelNumber = analogChannelsStored(analogChannelPos);
            else
                analogChannelNumber = 0;
            end
            channelNumber = [analogChannelNumber spikeChannelNumber];
            
        elseif analysisType == 2 || analysisType == 3 % Do once Spike data is available
            channelPos = get(hNeuralChannel,'val');
            channelNumber = neuralChannelsStored(channelPos);
            unitID = SourceUnitIDs(channelPos);
            plotSpikeData1Channel(plotHandles,channelNumber,s,f,o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,p);
            plotSpikeData1Parameter1Channel(hTemporalFreqPlot,channelNumber,a,e,s,f,o,c,[],folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hContrastPlot,channelNumber,a,e,s,f,o,[],t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hOrientationPlot,channelNumber,a,e,s,f,[],c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hSpatialFreqPlot,channelNumber,a,e,s,[],o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
            plotSpikeData1Parameter1Channel(hSigmaPlot,channelNumber,a,e,[],f,o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
        else
            analogChannelPos = get(hAnalogChannel,'val');
            rfMapVals = plotLFPData1Channel(plotHandles,analogChannelPos,s,f,o,...
                analysisType,lfpData,plotColor,blRange,stRange,p);
            plotLFPData1Parameter1Channel(hOrientationPlot,analogChannelPos,a,e,s,f,[],...
                analysisType,lfpData,plotColor,blRange,stRange,p);
            plotLFPData1Parameter1Channel(hSpatialFreqPlot,analogChannelPos,a,e,s,[],o,...
                analysisType,lfpData,plotColor,blRange,stRange,p);
            plotLFPData1Parameter1Channel(hSigmaPlot,analogChannelPos,a,e,[],f,o,...
                analysisType,lfpData,plotColor,blRange,stRange,p);

            channelNumber = analogChannelPos;
%             if analogChannelPos<=length(analogChannelsStored)
%                 channelNumber = analogChannelsStored(analogChannelPos);
%             else
%                 channelNumber = 0;
%             end
            
%             if ~isempty(rfMapVals)
%                 if (length(p.aValsUnique)==1) || (length(p.eValsUnique)==1)
%                     disp('Not enough data to plot RF center...')
%                 else
%                     plotRFMaps(hRFMapPlot,hcenterRFMapPlot,rfMapVals,p.aValsUnique,p.eValsUnique,plotColor,holdOnState);
%                 end
%             end
        end

        if analysisType<=3  % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType <=5
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        else
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        end

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(hOrientationPlot,xMin,xMax,getYLims(hOrientationPlot));
        rescaleData(hSpatialFreqPlot,xMin,xMax,getYLims(hSpatialFreqPlot));
        rescaleData(hSigmaPlot,xMin,xMax,getYLims(hSigmaPlot));
        showElectrodeLocationsMayo([],channelNumber,plotColor,hElectrodes,holdOnState,0);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        
        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        yLims = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,xMin,xMax,yLims);
        rescaleData(hOrientationPlot,xMin,xMax,yLims);
        rescaleData(hSpatialFreqPlot,xMin,xMax,yLims);
        rescaleData(hSigmaPlot,xMin,xMax,yLims);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType==6
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        else    
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(hOrientationPlot,xMin,xMax,getYLims(hOrientationPlot));
        rescaleData(hSpatialFreqPlot,xMin,xMax,getYLims(hSpatialFreqPlot));
        rescaleData(hSigmaPlot,xMin,xMax,getYLims(hSigmaPlot));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
        
        holdOnGivenPlotHandle(plotHandles,holdOnState);
        holdOnGivenPlotHandle(hOrientationPlot,holdOnState);
        holdOnGivenPlotHandle(hSpatialFreqPlot,holdOnState);
        holdOnGivenPlotHandle(hSigmaPlot,holdOnState);
        
        if holdOnState
            set(hElectrodes,'Nextplot','add');
        else
            set(hElectrodes,'Nextplot','replace');
        end

        function holdOnGivenPlotHandle(plotHandles,holdOnState)
            
            [numRows,numCols] = size(plotHandles);
            if holdOnState
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','add');

                    end
                end
            else
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','replace');
                    end
                end
            end
        end 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        
        claGivenPlotHandle(plotHandles);
        claGivenPlotHandle(hOrientationPlot);
        claGivenPlotHandle(hSpatialFreqPlot);
        claGivenPlotHandle(hSigmaPlot);
        
        cla(hRFMapPlot);cla(hcenterRFMapPlot);
        
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function rfMapVals = plotLFPData1Channel(plotHandles,channelNumber,s,f,o,...
analysisType,lfpData,plotColor,blRange,stRange,p)

titleFontSize = 10;

[numRows,numCols] = size(plotHandles);

% Get bad trials
% badTrialFile = fullfile(folderSegment,'badTrials.mat');
% if ~exist(badTrialFile,'file')
%     disp('Bad trial file does not exist...');
%     badTrials=[];
% else
%     badTrials = loadBadTrials(badTrialFile);
%     disp([num2str(length(badTrials)) ' bad trials']);
% end
timeVals=lfpData.timeVals;
analogData=squeeze(lfpData.segmentedLFPData(channelNumber,:,:));
badTrials=[];
rfMapVals = zeros(numRows,numCols);
for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = p.parameterCombinations{a,e,s,f,o};
        goodPos = setdiff(goodPos,badTrials);
      
        if isempty(goodPos)
            disp('No entries for this combination..');
        else
            disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);
    
            Fs = round(1/(timeVals(2)-timeVals(1)));
            if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
                disp('baseline and stimulus ranges are not the same');
            else
                range = blRange;
                rangePos = round(diff(range)*Fs);
                blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
                stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
                xs = 0:1/diff(range):Fs-1/diff(range);
            end

            if analysisType == 1        % compute ERP
                clear erp
                erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
                plot(plotHandles(i,j),timeVals,erp,'color',plotColor);
                
                rfMapVals(e,a) = rms(erp(stPos));

            elseif analysisType == 2  ||   analysisType == 3 % compute Firing rates
                disp('Use plotSpikeData instead of plotLFPData...');
            else
                fftBL = abs(fft(analogData(goodPos,blPos),[],2));
                fftST = abs(fft(analogData(goodPos,stPos),[],2));

                if analysisType == 4
                    plot(plotHandles(i,j),xs,log10(mean(fftBL)),'g');
                    set(plotHandles(i,j),'Nextplot','add');
                    plot(plotHandles(i,j),xs,log10(mean(fftST)),'k');
                    set(plotHandles(i,j),'Nextplot','replace');
                end

                if analysisType == 5
                    plot(plotHandles(i,j),xs,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                end
            end
            
            % Display title
            if (i==1)
                if (j==1)
                    title(plotHandles(i,j),['Azi: ' num2str(p.aValsUnique(a))],'FontSize',titleFontSize);
                else
                    title(plotHandles(i,j),num2str(p.aValsUnique(a)),'FontSize',titleFontSize);
                end
            end
                
            if (j==numCols)
                if (i==1)
                 title(plotHandles(i,j),[{'Ele'} {num2str(p.eValsUnique(e))}],'FontSize',titleFontSize,...
                     'Units','Normalized','Position',[1.25 0.5]);
                else
                    title(plotHandles(i,j),num2str(p.eValsUnique(e)),'FontSize',titleFontSize,...
                     'Units','Normalized','Position',[1.25 0.5]);
                end
            end
        end
    end
end

if analysisType~=1
    rfMapVals=[];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLFPData1Parameter1Channel(plotHandles,channelNumber,a,e,s,f,o,...
analysisType,lfpData,plotColor,blRange,stRange,p)

titleFontSize = 10;

% For orientation selectivity analysis
timeForComputation = [40 100]/1000; % s
freqForComputation = [40 60]; % Hz

parameterCombinations = p.parameterCombinations;
aValsUnique = p.aValsUnique;
eValsUnique = p.eValsUnique;
sValsUnique = p.sValsUnique;
fValsUnique = p.fValsUnique;
oValsUnique = p.oValsUnique;

[~,numCols] = size(plotHandles);

% Get the data
timeVals=lfpData.timeVals;
analogData=squeeze(lfpData.segmentedLFPData(channelNumber,:,:));
badTrials=[];

% Get bad trials
% badTrialFile = fullfile(folderSegment,'badTrials.mat');
% if ~exist(badTrialFile,'file')
%     disp('Bad trial file does not exist...');
%     badTrials=[];
% else
%     badTrials = loadBadTrials(badTrialFile);
%     disp([num2str(length(badTrials)) ' bad trials']);
% end

% Out of a,e,s,f,o,c and t only one parameter is empty
if isempty(a)
    aList = 1:length(aValsUnique); 
    eList = e+zeros(1,numCols); 
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    titleParam = 'Azi: ';
    titleList = aValsUnique;
end

if isempty(e) 
    aList = a+zeros(1,numCols);
    eList = 1:length(eValsUnique);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    titleParam = 'Ele: ';
    titleList = eValsUnique;
end

if isempty(s) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = 1:length(sValsUnique);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    titleParam = 'Sigma: ';
    titleList = sValsUnique;
end

if isempty(f)
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = 1:length(fValsUnique); 
    oList = o+zeros(1,numCols);
    titleParam = 'SF: ';
    titleList = fValsUnique;
end

if isempty(o) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = 1:length(oValsUnique);
    titleParam = 'Ori: ';
    titleList = oValsUnique;
end

% Main loop
computationVals=zeros(1,numCols);
for j=1:numCols
    clear goodPos
    goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j)};
    goodPos = setdiff(goodPos,badTrials);

    if isempty(goodPos)
        disp('No entries for this combination..')
    else
        disp(['pos=' num2str(j) ',n=' num2str(length(goodPos))]);

        Fs = round(1/(timeVals(2)-timeVals(1)));
        if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
            disp('baseline and stimulus ranges are not the same');
        else
            range = blRange;
            rangePos = round(diff(range)*Fs);
            blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
            stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
            xs = 0:1/diff(range):Fs-1/diff(range);
        end
        
        xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
        freqComputation = intersect(find(xs>=freqForComputation(1)),find(xs<=freqForComputation(2)));

        if analysisType == 1        % compute ERP
            clear erp
            erp = mean(analogData(goodPos,:),1);
            plot(plotHandles(j),timeVals,erp,'color',plotColor);
            
            if isempty(o) % Orientation tuning
                computationVals(j) = abs(min(erp(xsComputation)));
            end

        elseif analysisType == 2 || analysisType == 3   % compute Firing rates
            disp('Use plotSpikeData instead of plotLFPData...');
        else
            fftBL = abs(fft(analogData(goodPos,blPos),[],2));
            fftST = abs(fft(analogData(goodPos,stPos),[],2));

            if analysisType == 4
                plot(plotHandles(j),xs,log10(mean(fftBL)),'g');
                set(plotHandles(j),'Nextplot','add');
                plot(plotHandles(j),xs,log10(mean(fftST)),'k');
                set(plotHandles(j),'Nextplot','replace');
            end

            if analysisType == 5
                plot(plotHandles(j),xs,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
            end
            
            if isempty(o) % Orientation tuning
                computationVals(j) = max(mean(fftST(:,freqComputation),1));
            end
        end

        % Display title
        if (j==1)
            title(plotHandles(j),[titleParam num2str(titleList(j))],'FontSize',titleFontSize);
        else
            title(plotHandles(j),num2str(titleList(j)),'FontSize',titleFontSize);
        end
    end
end

% Orientation tuning
if isempty(o)
    if analysisType == 1
        disp(['Orientation selectivity values calculated between ' num2str(timeForComputation(1)) '-' num2str(timeForComputation(2))  ' s']);
    else
        disp(['Orientation selectivity values calculated between ' num2str(freqForComputation(1)) '-' num2str(freqForComputation(2))  ' Hz']);
    end
    disp(['orientation values: ' num2str(computationVals)]);
    [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique);
    disp(['prefOri: ' num2str(round(prefOrientation)) ', sel: ' num2str(orientationSelectivity)]);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Channel(plotHandles,channelNumber,s,f,o,c,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,folderName)
titleFontSize = 12;

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the data
clear signal spikeData
load(fullfile(folderSpikes,['elec' num2str(channelNumber) '_SID' num2str(unitID) '.mat']));

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s,f,o,c,t};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);
            
            if analysisType == 2
                [psthVals,xs] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
                plot(plotHandles(i,j),xs,psthVals,'color',plotColor);
            else
                X = spikeData(goodPos);
                axes(plotHandles(i,j)); %#ok<LAXES>
                rasterplot(X,1:length(X),plotColor);
            end
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Parameter1Channel(plotHandles,channelNumber,a,e,s,f,o,c,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,folderName)
titleFontSize = 12;

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);
[~,numCols] = size(plotHandles);

% Get the data
clear signal spikeData
load(fullfile(folderSpikes,['elec' num2str(channelNumber) '_SID' num2str(unitID) '.mat']));

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

% Out of a,e,s,f,o,c ant t, only one parameter is empty
if isempty(a)
    aList = 1:length(aValsUnique); 
    eList = e+zeros(1,numCols); 
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Azi: ';
    titleList = aValsUnique;
end

if isempty(e) 
    aList = a+zeros(1,numCols);
    eList = 1:length(eValsUnique);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Ele: ';
    titleList = eValsUnique;
end

if isempty(s) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = 1:length(sValsUnique);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Sigma: ';
    titleList = sValsUnique;
end

if isempty(f)
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = 1:length(fValsUnique); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'SF: ';
    titleList = fValsUnique;
end

if isempty(o) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = 1:length(oValsUnique);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Ori: ';
    titleList = oValsUnique;
end

if isempty(c) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = 1:length(cValsUnique);
    tList = t+zeros(1,numCols);
    titleParam = 'Con: ';
    titleList = cValsUnique;
end

if isempty(t) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols);
    tList = 1:length(tValsUnique);
    titleParam = 'TF: ';
    titleList = tValsUnique;
end

% Plot

for j=1:numCols
    %a = j;
    clear goodPos
    goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
    goodPos = setdiff(goodPos,badTrials);

    if isempty(goodPos)
        disp('No entries for this combination..')
    else
        disp(['pos=' num2str(j) ',n=' num2str(length(goodPos))]);
        if analysisType == 2
            [psthVals,xs] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
            plot(plotHandles(j),xs,psthVals,'color',plotColor);
        else
            X = spikeData(goodPos);
            axes(plotHandles(j)); %#ok<LAXES>
            rasterplot(X,1:length(X),plotColor);
        end
    end

    % Display title
    if (j==1)
        title(plotHandles(j),[titleParam num2str(titleList(j))],'FontSize',titleFontSize);
    else
        title(plotHandles(j),num2str(titleList(j)),'FontSize',titleFontSize);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,blRange,stRange,folderName,staLen,removeMeanSTA)

titleFontSize = 12;

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the analog data
clear signal analogData
load(fullfile(folderLFP,analogChannelString));

% Get the spike data
clear signal spikeData
load(fullfile(folderSpikes,['elec' num2str(spikeChannelNumber) '_SID' num2str(unitID) '.mat']));

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

staTimeLims{1} = blRange;
staTimeLims{2} = stRange;

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s,f,o,c,t};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            goodSpikeData = spikeData(goodPos);
            goodAnalogSignal = analogData(goodPos,:);
            [staVals,numberOfSpikes,xsSTA] = getSTA(goodSpikeData,goodAnalogSignal,staTimeLims,timeVals,staLen,removeMeanSTA);
            
            disp([num2str(i) ' ' num2str(j) ', numStim: ' num2str(length(goodPos)) ', numSpikes: ' num2str(numberOfSpikes)]);
            if ~isempty(staVals{1})
                plot(plotHandles(i,j),xsSTA,staVals{1},'color',plotColors{1});
            end
            set(plotHandles(i,j),'Nextplot','add');
            if ~isempty(staVals{2})
                plot(plotHandles(i,j),xsSTA,staVals{2},'color',plotColors{2});
            end
            set(plotHandles(i,j),'Nextplot','replace');
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotSTA1Parameter1Channel(hOrientationPlot,analogChannelNumber,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
%                 a,e,s,f,o,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRFMaps(hRFMapPlot,hcenterRFMapPlot,rfMapVals,aValsUnique,eValsUnique,plotColor,holdOnState)

plotSimple = 0;% Just compute the mean and Variance

if plotSimple
    [aziCenter,eleCenter] = getRFcenterSimple(aValsUnique,eValsUnique,rfMapVals); %#ok<*UNRCH>
else
    outParams = getRFcenter(aValsUnique,eValsUnique,rfMapVals);
    aziCenter = outParams(1); eleCenter = outParams(2);
    RFSize = sqrt((outParams(3)^2+outParams(4)^2)/2);
end

if plotSimple
    set(hRFMapPlot,'visible','off')
else
    % Plot the gaussian
    dX = (aValsUnique(end)-aValsUnique(1))/100;
    dY = (eValsUnique(end)-eValsUnique(1))/100;
    [~,outVals,boundaryX,boundaryY] = gauss2D(outParams,aValsUnique(1):dX:aValsUnique(end),eValsUnique(1):dY:eValsUnique(end));
    pcolor(hRFMapPlot,aValsUnique(1):dX:aValsUnique(end),eValsUnique(1):dY:eValsUnique(end),outVals);
    shading(hRFMapPlot,'interp'); 
    set(hRFMapPlot,'Nextplot','add');
    plot(hRFMapPlot,boundaryX,boundaryY,'k');
    set(hRFMapPlot,'Nextplot','replace');
end

% Plot the center only
if holdOnState
    set(hcenterRFMapPlot,'Nextplot','add');
else
    set(hcenterRFMapPlot,'Nextplot','replace');
end
plot(hcenterRFMapPlot,aziCenter,eleCenter,[plotColor '+']);
title(hcenterRFMapPlot,['Loc: (' num2str(aziCenter) ',' num2str(eleCenter) '), Size: ' num2str(RFSize)]);
disp(['Loc: (' num2str(aziCenter) ',' num2str(eleCenter) '), Size: ' num2str(RFSize)]);
axis(hcenterRFMapPlot,[aValsUnique(1) aValsUnique(end) eValsUnique(1) eValsUnique(end)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function rescaleData(plotHandles,xMin,xMax,yLims)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        else
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = getStringFromValues(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end
function [outString,outArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums)
outString='';
count=1;
for i=1:length(analogChannelsStored)
    outArray{count} = ['elec' num2str(analogChannelsStored(i))]; %#ok<AGROW>
    outString = cat(2,outString,[outArray{count} '|']);
    count=count+1;
end
if ~isempty(analogInputNums)
    for i=1:length(analogInputNums)
        outArray{count} = ['ainp' num2str(analogInputNums(i))];
        outString = cat(2,outString,[outArray{count} '|']);
        count=count+1;
    end
end
end
function outString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs)
outString='';
for i=1:length(neuralChannelsStored)
    outString = cat(2,outString,[num2str(neuralChannelsStored(i)) ', SID ' num2str(SourceUnitIDs(i)) '|']);
end 
end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data

function badTrials = loadBadTrials(badTrialFile) %#ok<*STOUT>
load(badTrialFile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [aziCenter,eleCenter] = getRFcenterSimple(aValsUnique,eValsUnique,rfMapValues)

rfMapValues = rfMapValues-mean(mean(rfMapValues));
% Compute the mean and variance
[maxEnonUnique,maxAnonUnique] = (find(rfMapValues==max(max(rfMapValues))));
maxE = maxEnonUnique(ceil(length(maxEnonUnique)/2));
maxA = maxAnonUnique(ceil(length(maxAnonUnique)/2));

%[maxE maxA]
aziCenter = sum(rfMapValues(maxE,:).*aValsUnique)/sum(rfMapValues(maxE,:));
eleCenter = sum(rfMapValues(:,maxA)'.*eValsUnique)/sum(rfMapValues(:,maxA));
end
function [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique)
num=0;
den=0;

for j=1:length(oValsUnique)
    num = num+computationVals(j)*sind(2*oValsUnique(j));
    den = den+computationVals(j)*cosd(2*oValsUnique(j));
end

prefOrientation = 90*atan2(num,den)/pi;
orientationSelectivity = abs(den+1i*num)/sum(computationVals);

if prefOrientation<0
    prefOrientation = prefOrientation+180;
end
end