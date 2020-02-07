% Display All Channels
function displayAllChannelsMayo(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end

% Folders and file names
folderName = fullfile(folderSourceString,'Data','segmentedData',fileNameString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.175; panelStartHeight = 0.775; backgroundColor = 'w'; % Panels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Parameters panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hParameterPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.35 panelStartHeight 0.2 panelHeight]);

% Trial Outcome
paramsHeight=1/6;

% Trial Type
trialTypeString = [{'Normal'} {'Instruction'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-1*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Trial Type','FontSize',fontSizeSmall);
hTrialType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-1*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',trialTypeString,'FontSize',fontSizeSmall);

trialOutcomeString = [{'Hit'} {'Missed'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-2*paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Trial Outcome','FontSize',fontSizeSmall);
hTrialOutcome = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-2*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',trialOutcomeString,'FontSize',fontSizeSmall);

% CueType
cueTypeString = [{'Valid'} {'Invalid'} {'Neutral'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-3*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Cue Type','FontSize',fontSizeSmall);
hCueType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [0.5 1-3*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',cueTypeString,'FontSize',fontSizeSmall);

% AttendLoc
attendLocString = [{'0'} {'1'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-4*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','AttendLoc','FontSize',fontSizeSmall);
hAttendLocType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [0.5 1-4*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',attendLocString,'FontSize',fontSizeSmall);

% Analysis Interval
analysisIntervalString = [{'StimOnset'} {'TargetOnset'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-5*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);
hTimePeriodType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-5*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',analysisIntervalString ,'FontSize',fontSizeSmall);

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|FFT|delta FFT';
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-6*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-6*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 1/6; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.55 panelStartHeight 0.2 panelHeight]);

signalRange = [-0.25 0.5]; fftRange = [0 100];
baseline = [-0.25 0]; stimPeriod = [0.25 0.5];

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
    'Position',[0 1-2*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-2*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-2*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-4*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.15;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.75 panelStartHeight 0.2 panelHeight]);

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
% Show electrode array and bad channels
electrodeGridPos = [0.05 panelStartHeight 0.25 panelHeight];

% Bad channels List
badChannels=[];
hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,badChannels,'r',[],0,0);

% Get main plot and message handles
[~,~,electrodeArray] = electrodePositionOnGridMayo(1);
[numRows,numCols] = size(electrodeArray);
plotHandles = getPlotHandles(numRows,numCols,[0.05 0.05 0.9 0.7]);

hMessage = uicontrol('Unit','Normalized','Position',[0 0.95 1 0.05],...
    'Style','text','String',fileNameString,'FontSize',fontSizeLarge);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        tt=get(hTrialType,'val'); ttStr = trialTypeString{tt};
        t=get(hTrialOutcome,'val'); tStr = trialOutcomeString{t}(1);
        c=get(hCueType,'val'); cStr = cueTypeString{c}(1);
        a=get(hAttendLocType,'val'); aStr = attendLocString{a};
        tp=get(hTimePeriodType,'val'); tpStr = analysisIntervalString{tp};
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));

        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];

        % Load Data
        if strcmp(ttStr,'Instruction')
            appendStr = '_Instruct';
        else
            appendStr = '';
        end
        if strcmp(cStr,'N')
            fileNameSaveStringLFP = fullfile(folderName,[fileNameString tStr cStr '_' tpStr '_LFP' appendStr]);
            fileNameSaveStringSpikes = fullfile(folderName,[fileNameString tStr cStr '_' tpStr '_Spikes' appendStr]);
        else
            fileNameSaveStringLFP = fullfile(folderName,[fileNameString tStr aStr cStr '_' tpStr '_LFP' appendStr]);
            fileNameSaveStringSpikes = fullfile(folderName,[fileNameString tStr aStr cStr '_' tpStr '_Spikes' appendStr]);
        end
        
        lfpData=load(fileNameSaveStringLFP);
        spikeData=load(fileNameSaveStringSpikes);
        
        set(hMessage,'String',[num2str(size(lfpData.segmentedLFPData,2)) ' stimuli found' ]);
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end
        
        if analysisType == 2 % Spikes
            plotSpikeData(plotHandles,spikeData,plotColor,xRange);
        else
            plotLFPData(plotHandles,lfpData,analysisType,plotColor,blRange,stRange);
        end

        yRange = getYLims(plotHandles);
        set(hYMin,'String',num2str(yRange(1))); set(hYMax,'String',num2str(yRange(2)));
        rescaleData(plotHandles,[xRange yRange]);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end

        yRange = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,[xRange yRange]);

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end

        yRange = getYLims(plotHandles);
        rescaleData(plotHandles,[xRange yRange]);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');

        [numRow,numCol] = size(plotHandles);

        if holdOnState
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','add');
                end
            end
        else
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','replace');
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        [numRow,numCol] = size(plotHandles);
        for i=1:numRow
            for j=1:numCol
                cla(plotHandles(i,j));
            end
        end

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function plotLFPData(plotHandles,lfpData,analysisType,plotColor,blRange,stRange)

timeVals = lfpData.timeVals;

if analysisType>2 % FFT
    Fs = round(1/(timeVals(2)-timeVals(1)));
    blPos = find(timeVals>=blRange(1),1)+ (0:diff(blRange)*Fs-1);
    stPos = find(timeVals>=stRange(1),1)+ (0:diff(stRange)*Fs-1);
    
    xsBL = 0:1/(diff(blRange)):Fs-1/(diff(blRange));
    xsST = 0:1/(diff(stRange)):Fs-1/(diff(stRange));
end

numElectrodes = size(lfpData.segmentedLFPData,1);
for i=1:numElectrodes
    
    disp(['Plotting electrode ' num2str(i)]);
    
    % get position
    [row,column] = electrodePositionOnGridMayo(i);
    
    if analysisType == 1        % compute ERP
        erp = squeeze(mean(lfpData.segmentedLFPData(i,:,:),2)); 
        plot(plotHandles(row,column),timeVals,erp,'color',plotColor);
        
    elseif analysisType == 2    % compute Firing rates
        disp('Use plotSpikeData instead of plotLFPData...');
    else
        fftBL = squeeze(abs(fft(lfpData.segmentedLFPData(i,:,blPos),[],3)));
        fftST = squeeze(abs(fft(lfpData.segmentedLFPData(i,:,stPos),[],3)));
        
        if analysisType == 3
            plot(plotHandles(row,column),xsBL,log10(mean(fftBL)),'k--');
            hold(plotHandles(row,column),'on');
            plot(plotHandles(row,column),xsST,log10(mean(fftST)),plotColor);
            hold(plotHandles(row,column),'off');
        end
        
        if analysisType == 4
            if isequal(xsBL,xsST)
                plot(plotHandles(row,column),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                hold(plotHandles(row,column),'on');
                plot(plotHandles(row,column),xsBL,zeros(1,length(xsBL)),'k--');
                hold(plotHandles(row,column),'off');
            else
                disp('Choose same baseline and stimulus periods..');
            end
        end
    end
end
end
function plotSpikeData(plotHandles,spikeData,plotColor,xRange)

binWidthMS = 10;

numElectrodes = size(spikeData.segmentedSpikeData,1);

for i=1:numElectrodes
    [row,column] = electrodePositionOnGridMayo(i); % get position
    [psthVals,xs] = getPSTH(spikeData.segmentedSpikeData(i,:),binWidthMS,[xRange(1) xRange(2)]);
    plot(plotHandles(row,column),xs,psthVals,'color',plotColor);
end
end   

function yRange = getYLims(plotHandles,channelsStored)

if ~exist('channelsStored','var');    channelsStored=1:96;              end

% Initialize
yMin = inf; yMax = -inf;

for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    % get position
    [row,column] = electrodePositionOnGridMayo(channelNum);
    
    axis(plotHandles(row,column),'tight');
    tmpAxisVals = axis(plotHandles(row,column));
    if tmpAxisVals(3) < yMin
        yMin = tmpAxisVals(3);
    end
    if tmpAxisVals(4) > yMax
        yMax = tmpAxisVals(4);
    end
end
yRange = [yMin yMax];
end
function rescaleData(plotHandles,axisLims,channelsStored)

if ~exist('channelsStored','var');    channelsStored=1:96;              end

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    
    % get position
   [row,column] = electrodePositionOnGridMayo(channelNum);
    
    axis(plotHandles(row,column),axisLims);
    if (row==numRows && rem(column,2)==rem(numCols+1,2))
        if column~=1
            set(plotHandles(row,column),'YTickLabel',[],'fontSize',labelSize);
        end
    elseif (rem(row,2)==0 && column==1)
        set(plotHandles(row,column),'XTickLabel',[],'fontSize',labelSize);
    else
        set(plotHandles(row,column),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
    end
end

% Remove Labels on the four corners
set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colorString, colorNames] = getColorString
colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';
end