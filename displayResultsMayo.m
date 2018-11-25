function displayResultsMayo(folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.8; backgroundColor = 'w'; % Panels

% Figure
close all;
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1])

% Electrode Grid will now be shown according to session chosen
% Show electrodes 
% electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
% hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hParameterPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.25 panelStartHeight 0.25 panelHeight]);
paramsHeight=1/6;

% FileNameString
[fileNameStringAll,fileNameStringListAll,fileNameStringListArray] = getFileNameStringList;
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Session(s)','FontSize',fontSizeSmall);
hSession = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',fileNameStringAll,'FontSize',fontSizeSmall);

% Type of neuron - single unit (SUA), multi-unit (MUA) or all
neuronTypeString = [{'All'} {'SUA'} {'MUA'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-2*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Neuron Type','FontSize',fontSizeSmall);
hNeuronType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [0.5 1-2*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',neuronTypeString,'FontSize',fontSizeSmall);

% Population Type - either take all electrodes for analysis, or only that
% are stimulated by the stimulus
populationString = [{'All'} {'Stimulated'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-3*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Population','FontSize',fontSizeSmall);
hPopulationType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [0.5 1-3*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',populationString,'FontSize',fontSizeSmall);

trialOutcomeString = [{'Hit'} {'Missed'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-4*paramsHeight 0.5 paramsHeight],...
    'Style','text','String','Trial Outcome','FontSize',fontSizeSmall);
hTrialOutcome = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-4*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',trialOutcomeString,'FontSize',fontSizeSmall);

% Orientation Change
orientationChangeString=[{'All'} {'Selected'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-5*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Orientation Change','FontSize',fontSizeSmall);
hOriChangeType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-5*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',orientationChangeString ,'FontSize',fontSizeSmall);

% Analysis Interval
analysisIntervalString = [{'Baseline'} {'StimOnset'} {'TargetOnset'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-6*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);
hTimePeriodType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-6*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',analysisIntervalString ,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','X and Y Limits','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.5 panelStartHeight 0.25 panelHeight]);
timingHeight = 1/6; 

signalRange = [-0.25 0.5]; fftRange = [0 100];

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

% Firing Range
frRange = [0 100];
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-4*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Firing Rate (s/s)','FontSize',fontSizeSmall);
hFRRangeMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(frRange(1)),'FontSize',fontSizeSmall);
hFRRangeMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(frRange(2)),'FontSize',fontSizeSmall);

% ERP Range
erpRange = [-100 100];
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','ERP Range (uV)','FontSize',fontSizeSmall);
hERPRangeMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(erpRange(1)),'FontSize',fontSizeSmall);
hERPRangeMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(erpRange(2)),'FontSize',fontSizeSmall);

% FFT Y Range
fftYRange = [2 5];
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','LogFFT (y-axis)','FontSize',fontSizeSmall);
hFFTYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftYRange(1)),'FontSize',fontSizeSmall);
hFFTYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftYRange(2)),'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.75 panelStartHeight 0.2 panelHeight]);
plotOptionsHeight = 1/6;

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 5*plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeSmall);
hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 5*plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleXY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots

hPSTH = getPlotHandles(3,1,[0.16 0.05 0.1 0.7],0,0.05,0); linkaxes(hPSTH);
hERP = getPlotHandles(3,1,[0.28 0.05 0.1 0.7],0,0.05,0); linkaxes(hERP);
hFFT = getPlotHandles(3,1,[0.4 0.05 0.1 0.7],0,0.05,0); linkaxes(hFFT);
hFFC = getPlotHandles(3,1,[0.52 0.05 0.1 0.7],0,0.05,0); linkaxes(hFFC);
hSFC = getPlotHandles(3,1,[0.64 0.05 0.1 0.7],0,0.05,0); linkaxes(hSFC);


hBarFR = getPlotHandles(3,1,[0.76 0.05 0.04 0.7],0,0.05,0); linkaxes(hBarFR);
hBarRsc = getPlotHandles(3,1,[0.82 0.05 0.04 0.7],0,0.05,0); linkaxes(hBarRsc);
hBarAlpha = getPlotHandles(3,1,[0.88 0.05 0.04 0.7],0,0.05,0); linkaxes(hBarAlpha);
hBarSSVEP = getPlotHandles(3,1,[0.94 0.05 0.04 0.7],0,0.05,0); linkaxes(hBarSSVEP);

hBehavior(1) = subplot('Position',[0.03 0.3 0.1 0.45]);
hBehavior(2) = subplot('Position',[0.03 0.05 0.1 0.2]); linkaxes(hBehavior);

colorNamesSides = 'cm';

% Plotting functions
    function plotData_Callback(~,~)
        
        s=get(hSession,'val'); fileNameStringTMP = fileNameStringListArray{s}; SessionIDString = fileNameStringListAll{s};
        n=get(hNeuronType,'val'); neuronType = neuronTypeString{n};
        p=get(hPopulationType,'val'); populationType = populationString{p};
        t=get(hTrialOutcome,'val'); tStr = trialOutcomeString{t}(1);
        o=get(hOriChangeType,'val'); oStr=orientationChangeString{o};
        tp=get(hTimePeriodType,'val'); tpStr = analysisIntervalString{tp};
        
        colorNamesAttCue = colorNames{get(hChooseColor,'val')};
        
        % Show electrodes
%         if length(fileNameStringTMP)<=12
        electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
        hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,SessionIDString{1});
%         else
%         electrodeGridPos1 = [0.05 panelStartHeight+panelHeight/2 0.2 panelHeight/2]; 
%         hElectrodes1 = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,fileNameStringTMP);
%         electrodeGridPos2 = [0.05 panelStartHeight-panelHeight/2 0.2 panelHeight/2];
%         hElectrodes2 = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,fileNameStringTMP);
%         end
        
        if strcmp(tpStr,'TargetOnset')
            tRange = [-0.5 0.1];
        else
            tRange = [-0.25 0.5];
        end
        
        % Set Taper for Coherency analysis by MT method
        tapers = [1 1];
        
        % Get data
        [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ffcData,ffPhiData,sfcData,sfPhiData,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getData(folderSourceString,fileNameStringTMP,neuronType,populationType,tStr,oStr,tpStr,tapers);
        
        for attCuePos=1:5
            plot(hBehavior(1),uniqueOrientationChangeDeg,perCorrect(attCuePos,:),'color',colorNamesAttCue(attCuePos,:),'marker','o'); axis tight;
            hold(hBehavior(1),'on');
        end

        for arraySide=1:2
            showElectrodeLocationsMayo([],electrodeArray{arraySide},colorNamesSides(arraySide),hElectrodes,1,0,SessionIDString{1}); % Show electrodes used for analysis
              
            for attCuePos=1:5
                plot(hPSTH(arraySide),xsFR,squeeze(mean(psthData{arraySide}(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos,:));
                hold(hPSTH(arraySide),'on');
                
                plot(hERP(arraySide),timeVals,squeeze(mean(erpData{arraySide}(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos,:));
                hold(hERP(arraySide),'on');
                
                plot(hFFT(arraySide),freqVals,squeeze(mean(fftData{arraySide}(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos,:));
                hold(hFFT(arraySide),'on');
            end
            makeBarPlot(hBarFR(arraySide),firingRates{arraySide},colorNamesAttCue);
            makeBarPlot(hBarAlpha(arraySide),alphaData{arraySide},colorNamesAttCue);
            makeBarPlot(hBarSSVEP(arraySide),ssvepData{arraySide},colorNamesAttCue);
        end
        
        % Rescale plots and set the x and y scales
        yLims = getYLims(hPSTH);
        axis(hPSTH(1),[tRange 0 yLims(2)]);
        set(hStimMin,'String',num2str(tRange(1))); set(hStimMax,'String',num2str(tRange(2))); % Set tRange 
        set(hFRRangeMax,'String',num2str(yLims(2)));
        
        yLims = getYLims(hERP); 
        axis(hERP(1),[tRange yLims]);
        set(hERPRangeMin,'String',num2str(yLims(1))); set(hERPRangeMax,'String',num2str(yLims(2)));
        
        yLims = getYLims(hFFT);
        axis(hFFT(1),[str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String')) yLims]);
        set(hFFTYMin,'String',num2str(yLims(1))); set(hFFTYMax,'String',num2str(yLims(2)));
        
        yLims = getYLims(hBarFR(1:2)); axis(hBarFR(1),[0 6 0 yLims(2)]);
        yLims = getYLims(hBarAlpha(1:2)); axis(hBarAlpha(1),[0 6 yLims]);
        yLims = getYLims(hBarSSVEP(1:2)); axis(hBarSSVEP(1),[0 6 yLims]);
        ylim(hBehavior(1),[0 1]);
        
        % Labels
        xlabel(hPSTH(3),'Time (s)'); title(hPSTH(1),'Firing Rate (spikes/s)');
        xlabel(hERP(3),'Time (s)'); title(hERP(1),'ERP (\muV)');
        xlabel(hFFT(3),'Frequency (Hz)'); title(hFFT(1),'Log FFT');
        xlabel(hFFC(3),'Frequency (Hz)'); title(hFFC(1),'FFC');
        xlabel(hSFC(3),'Frequency (Hz)'); title(hSFC(1),'SFC');
        
        title(hBarFR(1),'Firing Rate');
        title(hBarRsc(1),'r_s_c');
        title(hBarAlpha(1),'Alpha');
        title(hBarSSVEP(1),'SSVEP');
        
        title(hBehavior(1),'Behavior (Fraction Correct)');
        
        for arraySide=1:2
            ylabel(hPSTH(arraySide),['N=' num2str(length(electrodeArray{arraySide}))],'color',colorNamesSides(arraySide));
        end
        
        legend(hBehavior(1),'0V','1V','0I','1I','N','location','southeast');
    end
    function rescaleXY_Callback(~,~)
        tRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        fRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        
        axis(hPSTH(1),[tRange str2double(get(hFRRangeMin,'String')) str2double(get(hFRRangeMax,'String'))]);
        axis(hERP(1),[tRange str2double(get(hERPRangeMin,'String')) str2double(get(hERPRangeMax,'String'))]);
        axis(hFFT(1),[fRange str2double(get(hFFTYMin,'String')) str2double(get(hFFTYMax,'String'))]);
    end
    function holdOn_Callback(~,~)
        holdOnState = get(hHoldOn,'Value');
        
        holdOnGivenPlotHandle(hPSTH,holdOnState);
        holdOnGivenPlotHandle(hERP,holdOnState);
        holdOnGivenPlotHandle(hFFT,holdOnState);
        holdOnGivenPlotHandle(hBarFR,holdOnState);
        holdOnGivenPlotHandle(hBarAlpha,holdOnState);
        holdOnGivenPlotHandle(hBarSSVEP,holdOnState);
        holdOnGivenPlotHandle(hBehavior,holdOnState);
        
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
    function cla_Callback(~,~)
        
        claGivenPlotHandle(hPSTH);
        claGivenPlotHandle(hERP);
        claGivenPlotHandle(hFFT);
        claGivenPlotHandle(hBarFR);
        claGivenPlotHandle(hBarAlpha);
        claGivenPlotHandle(hBarSSVEP);
        claGivenPlotHandle(hBehavior);
        
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
function [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ffcData,ffPhiData,sfcData,sfPhiData,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getData(folderSourceString,fileNameStringTMP,neuronType,populationType,tStr,oStr,tpStr,tapers)

numDatasets = length(fileNameStringTMP);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ffcData,ffPhiData,sfcData,sfPhiData,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getDataSingleSession(folderSourceString,fileNameStringTMP{1},neuronType,populationType,tStr,oStr,tpStr,tapers); % First session

if length(fileNameStringTMP)>1
    for i=2:numDatasets
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [psthDataTMP,xsFRTMP,firingRatesTMP,erpDataTMP,timeValsTMP,fftDataTMP,freqValsTMP,alphaDataTMP,ssvepDataTMP,ffcDataTMP,ffPhiDataTMP,sfcDataTMP,sfPhiDataTMP,electrodeArrayTMP,perCorrectTMP,uniqueOrientationChangeDegTMP] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},neuronType,populationType,tStr,oStr,tpStr,tapers);
        
        perCorrect = perCorrect + perCorrectTMP;
        uniqueOrientationChangeDeg = uniqueOrientationChangeDeg + uniqueOrientationChangeDegTMP;
        
        for k=1:2 % for each array side
            if isequal(xsFR,xsFRTMP)
                psthData{k} = cat(2,psthData{k},psthDataTMP{k});
                firingRates{k} = cat(2,firingRates{k},firingRatesTMP{k});
            else
                error('xsFR does not match');
            end
            if isequal(timeVals,timeValsTMP)
                erpData{k} = cat(2,erpData{k},erpDataTMP{k});
            else
                error('timeVals do not match');
            end
            if isequal(freqVals,freqValsTMP)
                fftData{k} = cat(2,fftData{k},fftDataTMP{k});
                alphaData{k} = cat(2,alphaData{k},alphaDataTMP{k});
                ssvepData{k} = cat(2,ssvepData{k},ssvepDataTMP{k});
                ffcData{k} = cat(2,ffcData{k},ffcDataTMP{k});
                ffPhiData{k} = cat(2,ffPhiData{k},ffPhiDataTMP{k});
                sfcData{k} = cat(2,sfcData{k},sfcDataTMP{k});
                sfPhiData{k} = cat(2,sfPhiData{k},sfPhiDataTMP{k});
            else
                error('freqVals do not match');
            end
            
            
            electrodeArray{k} = cat(1,electrodeArray{k},electrodeArrayTMP{k});
        end
    end
    perCorrect = perCorrect/numDatasets;
    uniqueOrientationChangeDeg = uniqueOrientationChangeDeg/numDatasets;
end
end
function [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ffcData,ffPhiData,sfcData,sfPhiData,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getDataSingleSession(folderSourceString,fileNameString,neuronType,populationType,tStr,oStr,tpStr,tapers)

folderSave = fullfile(folderSourceString,'Data','savedData');
makeDirectory(folderSave);

fileToSave = fullfile(folderSave,[fileNameString neuronType populationType tStr oStr tpStr '.mat']);
coherencyFileToSave = fullfile(folderSave,[fileNameString neuronType populationType tStr oStr tpStr 'tapers_',tapers(1) '.mat']);

if exist(fileToSave,'file')&& exist(coherencyFileToSave,'file')
    disp(['Loading file ' fileToSave]);
    disp(['Loading file ' coherencyFileToSave]);
    load(fileToSave);load(coherencyFileToSave);
else
    
    binWidthMS = 10;
    alphaRangeHz = [8 12]; ssvepFreqHz = 20;
    
    folderName = fullfile(folderSourceString,'Data','segmentedData',fileNameString);
    % get Good Electrodes
    electrodeArray = getGoodElectrodes(folderSourceString,fileNameString,neuronType,populationType);
    % get Good Electrode Pairs
    [electrodepairsWithinHemisphere, electrodePairsAcrossHemispheres] = getGoodElectrodePairs(electrodeArray);

    
    % Get Spike and LFP data for 5 conditions - 0V, 1V, 0I, 1I and N
    attCueList = [{'0V'} {'1V'} {'0I'} {'1I'} {'N'}];
    numConditions = length(attCueList);
    
    for i=1:numConditions
        if strcmp(tpStr,'Baseline')
            timeRange = [-0.25 0];
            lfpData{i}=load(fullfile(folderName,[fileNameString tStr attCueList{i} '_StimOnset_LFP']));
            spikeData{i}=load(fullfile(folderName,[fileNameString tStr attCueList{i} '_StimOnset_Spikes']));
        elseif strcmp(tpStr,'StimOnset')
            timeRange = [0.25 0.5];
            lfpData{i}=load(fullfile(folderName,[fileNameString tStr attCueList{i} '_StimOnset_LFP']));
            spikeData{i}=load(fullfile(folderName,[fileNameString tStr attCueList{i} '_StimOnset_Spikes']));
        elseif strcmp(tpStr,'TargetOnset')
            timeRange = [-0.5 0];
            lfpData{i}=load(fullfile(folderName,[fileNameString tStr attCueList{i} '_TargetOnset_LFP']));
            spikeData{i}=load(fullfile(folderName,[fileNameString tStr attCueList{i} '_TargetOnset_Spikes']));
        end
    end
    
    [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString);
    if strcmp(tStr,'H')
        goodList = [goodIndexList(1:4) goodIndexList(9)]; % Hit conditions for 0V, 1V, 0I, 1I and M
    else
        goodList = [goodIndexList(5:8) goodIndexList(10)]; % Miss conditions for 0V, 1V, 0I, 1I and M
    end
    
    for i=1:numConditions
        if strcmp(oStr,'Selected') % only select trials for which 2nd and 3rd orientation changes occured
            allOrientationsThisCondition = orientationChangeDeg(goodList{i});
            selectedPos = sort([find(allOrientationsThisCondition==uniqueOrientationChangeDeg(2)) find(allOrientationsThisCondition==uniqueOrientationChangeDeg(3))]);
            lfpData{i}.segmentedLFPData = lfpData{i}.segmentedLFPData(:,selectedPos,:);
            spikeData{i}.segmentedSpikeData = spikeData{i}.segmentedSpikeData(:,selectedPos);
        end
        disp([fileNameString tStr attCueList{i} ', Stimulus repeats: ' num2str(size(lfpData{i}.segmentedLFPData,2))]);
    end
    % Frequency analysis
    timeVals=lfpData{1}.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    pos = find(timeVals>=timeRange(1),1)+ (0:diff(timeRange)*Fs-1);
    freqVals = 0:1/(diff(timeRange)):Fs-1/(diff(timeRange));
       
    alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
    ssvepPos = find(freqVals==ssvepFreqHz);
    
%     tapers = [1 1];
    for i=1:2 % Each array side
        eList = electrodeArray{i};
        ePairList = electrodepairsWithinHemisphere{i};
        clear psthDataTMP firingRatesTMP erpDataTMP fftDataTMP alphaDataTMP ssvepDataTMP ffcTMP ffPhiTMP sfcTMP sfPhiTMP
        for j=1:numConditions
            % Firing Rates
            for k=1:length(eList)
                [psthDataTMP(j,k,:),xsFR] = getPSTH(spikeData{j}.segmentedSpikeData(eList(k),:),binWidthMS,[timeVals(1) timeVals(end)]);
                firingRatesTMP(j,k) = mean(getSpikeCounts(spikeData{j}.segmentedSpikeData(eList(k),:),timeRange))/diff(timeRange);
            end
            
            % ERP
            erpDataTMP(j,:,:) = squeeze(mean(lfpData{j}.segmentedLFPData(eList,:,:),2));
            
            % FFT
            fftDataTMP(j,:,:) = log10(squeeze(mean(abs(fft(lfpData{j}.segmentedLFPData(eList,:,pos),[],3)),2)));
            alphaDataTMP(j,:) = sum(fftDataTMP(j,:,alphaPos),3);
            ssvepDataTMP(j,:) = fftDataTMP(j,:,ssvepPos); %#ok<FNDSB>
            
            % Coherency and Spike-LFP Phase Analysis (FFC,SFC,sfPhi)
            disp(['numArray:' num2str(i) ', numCondition:' num2str(j)]);
            [ffcTMP(j,:,:,:),ffPhiTMP(j,:,:,:),freqFFC,sfcTMP(j,:,:,:),sfPhiTMP(j,:,:,:),freqSFC] = getCoherencyMeasures(lfpData{j}.segmentedLFPData(:,:,pos),spikeData{j}.segmentedSpikeData,ePairList,tapers,lfpData{j}.timeVals(pos),timeRange,freqVals);
            
        end
        
        psthData{i} = psthDataTMP;
        firingRates{i} = firingRatesTMP;
        erpData{i} = erpDataTMP;
        fftData{i} = fftDataTMP;
        alphaData{i} = alphaDataTMP;
        ssvepData{i} = ssvepDataTMP;
        
        ffcData{i} = ffcTMP;
        ffPhiData{i} = ffPhiTMP;
        sfcData{i} = sfcTMP;
        sfPhiData{i} = sfPhiTMP;
    end
    
    % Save data
    save(fileToSave,'psthData','xsFR','firingRates','erpData','timeVals','fftData','freqVals','alphaData','ssvepData','electrodeArray','perCorrect','uniqueOrientationChangeDeg');
    save(coherencyFileToSave,'ffcData','ffPhiData','sfcData','sfPhiData','freqVals');
end
end
function [colorString, colorNames] = getColorString
colorNames{1} = [0 0 1; 1 0 0; 0 1 1; 1 0 1; 0 1 0]; colorString{1} = 'BlueRedCyanMagentaGreen';
colorNames{2} = gray(6); colorString{2} = 'gray';
colorNames{3} = copper(6); colorString{3} = 'copper';
colorNames{4} = jet(5); colorString{4} = 'jet';
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
function electrodeArray = getGoodElectrodes(folderSourceString,fileNameString,neuronType,populationType)

suaCutoff = 3;
% Get sorting rating
sortRatings=sortRating(fileNameString,folderSourceString);
sua=intersect(find(sortRatings>0),find(sortRatings<=suaCutoff));
mua=find(sortRatings>suaCutoff);
all=union(sua,mua,'sorted');

[~,~,electrodeArrayPos]=electrodePositionOnGridMayo(1,fileNameString);

spkData{1}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H1V_StimOnset_Spikes']));
spkData{2}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H0V_StimOnset_Spikes']));
 
 if strcmp(neuronType,'SUA')
    electrodeArray{1}=intersect(sua,electrodeArrayPos(:,8:13)); % Right Array
    electrodeArray{2}=intersect(sua,electrodeArrayPos(:,1:6)); % Left Array

 elseif strcmp(neuronType,'MUA')
    electrodeArray{1}=intersect(mua,electrodeArrayPos(:,8:13));
    electrodeArray{2}=intersect(mua,electrodeArrayPos(:,1:6));
 
 elseif strcmp(neuronType,'All')
     electrodeArray{1}=intersect(all,electrodeArrayPos(:,8:13));
     electrodeArray{2}=intersect(all,electrodeArrayPos(:,1:6));
 end
 
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
function [electrodepairsWithinHemisphere, electrodePairsAcrossHemispheres] = getGoodElectrodePairs(electrodeArray)

% get good Electrode pairs - within hemisphere
for j=1:2 % each electrode array
    electrodepairsWithinHemisphere{j} = combnk(electrodeArray{j},2); % Getting different electrode pairs for FFC and SFC using MATLAB combnk function
end

% get good Electrode pairs - across hemispheres
electrodePairsAcrossHemispheres = setdiff(combnk([electrodeArray{1}; electrodeArray{2}],2),[electrodepairsWithinHemisphere{1};electrodepairsWithinHemisphere{2}],'rows');
    
end
function [ffc,ffPhi,freqFFC,sfc,sfPhi,freqSFC,N] = getCoherencyMeasures(lfpData,spikeData,electrodePair,tapers,timeVals,timeRange,freqVals)

% Set up MT
Fs              = round(1/(timeVals(2)-timeVals(1)));
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 1;

% spk = convertSpikeTimes2Bins(spikeData,timeRange,1000/Fs);
N = size(lfpData,2);

disp('Working on FFC Data')
% Field-Field coherence
for i=1:size(electrodePair,1)
clear lfp1 lfp2
disp(['FFC ElectrodePair:',num2str(i)]);
lfp1 = squeeze(lfpData(electrodePair(i,1),:,:));
lfp2 = squeeze(lfpData(electrodePair(i,2),:,:));

[ffc(i,:,:),ffPhi(i,:,:),S12_ffc(i,:,:),S1_ffc(i,:,:),S2_ffc(i,:,:),freqFFC]=coherencyc(lfp1',lfp2',params); %#ok<ASGLU,*AGROW>
end

% if params.trialave = 0;
% [ppc_ffc,PLF_ffc,without_phaseCovariation_ffc,traditional_ffc]= computeCoherencyFromSpectrum(S1_ffc,S2_ffc,S12_ffc,params);
% ffc = traditional_ffc;
% end
    

disp('Working on SFC Data')
% Spike-Field coherence
for i=1:size(electrodePair,1)
clear lfp spk
disp(['SFC ElectrodePair:',num2str(i)]);
lfp = squeeze(lfpData(electrodePair(i,1),:,:));
spk = convertSpikeTimes2Bins(spikeData(electrodePair(i,2),:,:),timeRange,1000/Fs);
[sfc(i,:,:),sfPhi(i,:,:),S12_sfc(i,:,:),S1_sfc(i,:,:),S2_sfc(i,:,:),freqSFC]=coherencycpb(lfp',spk,params); %#ok<ASGLU,*AGROW>
end

% Sanity Check
if isequal(freqFFC,freqSFC,freqVals)
else
    error('freqVals from fft & multitaper do not match!')
end

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
function makeBarPlot(h,data,colorNames)

N = size(data,2);
mData = mean(data,2);
semData = std(data,[],2)/sqrt(N);

for i=1:size(data,1)
    plot(h,i,mData(i),'color',colorNames(i,:),'marker','o');
    hold(h,'on');
    errorbar(h,i,mData(i),semData(i),'color',colorNames(i,:));
end
set(h,'XTick',[],'XTicklabel',[]);
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