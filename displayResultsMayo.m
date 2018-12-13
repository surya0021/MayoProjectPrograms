function displayResultsMayo(folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end
close all;

% Figure 1
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 1

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.8; backgroundColor = 'w'; % Panels

% Electrode Grid will now be shown according to session chosen
% Show electrodes
electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,'blank');

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
analysisIntervalString = [{'Baseline'} {'StimOnset'} {'TargetOnset_250ms'} {'TargetOnset_500ms'}];
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

hShowAbsoluteVals = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','Show Absolute Measures','FontSize',fontSizeMedium);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1*plotOptionsHeight 0.5 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0.5 1*plotOptionsHeight 0.5 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 0.5 plotOptionsHeight], ...
    'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleXY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0.5 0 0.5 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots

% Figure 1
figure(1)
hBehavior(1) = subplot('Position',[0.05 0.3 0.2 0.45]);
hBehavior(2) = subplot('Position',[0.05 0.05 0.2 0.2]); linkaxes(hBehavior);

hPSTH = getPlotHandles(3,1,[0.3 0.05 0.1 0.7],0,0.05,0); linkaxes(hPSTH);
hERP = getPlotHandles(3,1,[0.45 0.05 0.1 0.7],0,0.05,0); linkaxes(hERP);
hFFT = getPlotHandles(3,1,[0.6 0.05 0.1 0.7],0,0.05,0); linkaxes(hFFT);

hBarFR = getPlotHandles(3,1,[0.75 0.05 0.05 0.7],0,0.05,0); linkaxes(hBarFR);
hBarAlpha = getPlotHandles(3,1,[0.82 0.05 0.05 0.7],0,0.05,0); linkaxes(hBarAlpha);
hBarSSVEP = getPlotHandles(3,1,[0.9 0.05 0.05 0.7],0,0.05,0); linkaxes(hBarSSVEP);

% Figure 2 (electrode pairwise analysis results)
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) %Figure 2
hFFC = getPlotHandles(4,1,[0.04 0.05 0.1 0.9],0,0.02,0); linkaxes(hFFC);
hSFC = getPlotHandles(4,1,[0.16 0.05 0.1 0.9],0,0.02,0); linkaxes(hSFC);
hAmpCorr = getPlotHandles(4,1,[0.28 0.05 0.1 0.9],0,0.02,0); linkaxes(hAmpCorr);
hffPPC = getPlotHandles(4,1,[0.4 0.05 0.1 0.9],0,0.02,0); linkaxes(hffPPC);
hsfPPC = getPlotHandles(4,1,[0.52 0.05 0.1 0.9],0,0.02,0); linkaxes(hsfPPC);

hBarRsc = getPlotHandles(4,1,[0.64 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarRsc);
hBarFFCAlpha = getPlotHandles(4,1,[0.69 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarFFCAlpha);
hBarFFCSSVEP = getPlotHandles(4,1,[0.74 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarFFCSSVEP);
hBarSFCAlpha = getPlotHandles(4,1,[0.79 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarSFCAlpha);
hBarSFCSSVEP = getPlotHandles(4,1,[0.84 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarSFCSSVEP);
hBarAmpCorrAlpha = getPlotHandles(4,1,[0.89 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarAmpCorrAlpha);
hBarAmpCorrSSVEP = getPlotHandles(4,1,[0.94 0.05 0.04 0.9],0,0.02,0); linkaxes(hBarAmpCorrSSVEP);

colorNamesSides = 'cmkk';

% Plotting functions
    function plotData_Callback(~,~)
        
        s=get(hSession,'val'); fileNameStringTMP = fileNameStringListArray{s}; SessionIDString = fileNameStringListAll{s};
        n=get(hNeuronType,'val'); neuronType = neuronTypeString{n};
        p=get(hPopulationType,'val'); populationType = populationString{p};
        t=get(hTrialOutcome,'val'); tStr = trialOutcomeString{t}(1);
        o=get(hOriChangeType,'val'); oStr=orientationChangeString{o};
        tp=get(hTimePeriodType,'val'); tpStr = analysisIntervalString{tp};
        
        colorNamesAttCue = colorNames{get(hChooseColor,'val')};
        showAbsoluteValsFlag = get(hShowAbsoluteVals,'val');
        
        % Show electrodes
        if strcmp(SessionIDString{1},'all (N=24)')
            hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,'blank');
        else
            electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
            hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0,SessionIDString{1});
        end
        
        if strcmp(tpStr(1:11),'TargetOnset')
            tRange = [-0.5 0.1];
        else
            tRange = [-0.25 0.5];
        end
        
        % Set Taper for Coherency analysis by MT method
        tapers = [1 1];
        
        % Get data
        [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ampCorrData,rSCData,ffcData,ffPhiData,sfcData,sfPhiData,freqValsMT,ffppcData,sfppcData,freqValsPPC,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getData(folderSourceString,fileNameStringTMP,neuronType,populationType,tStr,oStr,tpStr,tapers); %#ok<ASGLU>
        
        for attCuePos=1:5
            plot(hBehavior(1),uniqueOrientationChangeDeg,perCorrect(attCuePos,:),'color',colorNamesAttCue(attCuePos,:),'marker','o'); axis tight;
            hold(hBehavior(1),'on');
        end
        
        for arraySide=1:3
            if s==length(fileNameStringListAll)||arraySide==3
            else
                showElectrodeLocationsMayo([],electrodeArray{arraySide},colorNamesSides(arraySide),hElectrodes,1,0,SessionIDString{1}); % Show electrodes used for analysis
            end
            for attCuePos=1:5
                if attCuePos==3 || attCuePos==4
                    continue
                end
                
                plot(hPSTH(arraySide),xsFR,squeeze(mean(psthData{arraySide}(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos,:));
                hold(hPSTH(arraySide),'on');
                
                plot(hERP(arraySide),timeVals,squeeze(mean(erpData{arraySide}(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos,:));
                hold(hERP(arraySide),'on');
                
                clear dataFFT
                if showAbsoluteValsFlag
                    dataFFT = squeeze(fftData{arraySide}(attCuePos,:,:));
                else
                    dataFFT = squeeze(fftData{arraySide}(attCuePos,:,:)) - squeeze(fftData{arraySide}(5,:,:));
                end
                plotData(hFFT(arraySide),freqVals,dataFFT,colorNamesAttCue(attCuePos,:));
            end
            
            makeBarPlot(hBarFR(arraySide),firingRates{arraySide},colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarAlpha(arraySide),alphaData{arraySide},colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarSSVEP(arraySide),ssvepData{arraySide},colorNamesAttCue,showAbsoluteValsFlag);
        end
        
        alphaRangeHz = [8 12]; ssvepFreqHz = 20;
        alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
        ssvepPos = find(freqVals==ssvepFreqHz);
        
        for arraySide=1:4
            for attCuePos=1:5
                if attCuePos==3 || attCuePos==4
                    continue
                end
                clear dataFFC dataSFC dataAmpCorr dataffPPC datasfPPC
                if showAbsoluteValsFlag
                    dataFFC = squeeze(ffcData{arraySide}(attCuePos,:,:));
                    dataSFC = squeeze(sfcData{arraySide}(attCuePos,:,:));
                    dataAmpCorr = squeeze(ampCorrData{arraySide}(attCuePos,:,:));
                    dataffPPC = squeeze(ffppcData{arraySide}(attCuePos,:,:));
                    datasfPPC = squeeze(sfppcData{arraySide}(attCuePos,:,:));
                    
                else
                    dataFFC = squeeze(ffcData{arraySide}(attCuePos,:,:)) - squeeze(ffcData{arraySide}(5,:,:));
                    dataSFC = squeeze(sfcData{arraySide}(attCuePos,:,:)) - squeeze(sfcData{arraySide}(5,:,:));
                    dataAmpCorr = squeeze(ampCorrData{arraySide}(attCuePos,:,:)) - squeeze(ampCorrData{arraySide}(5,:,:));
                    dataffPPC = squeeze(ffppcData{arraySide}(attCuePos,:,:)) - squeeze(ffppcData{arraySide}(5,:,:));
                    datasfPPC = squeeze(sfppcData{arraySide}(attCuePos,:,:)) - squeeze(sfppcData{arraySide}(5,:,:));                    
                end
                plotData(hFFC(arraySide),freqValsMT,dataFFC,colorNamesAttCue(attCuePos,:));
                plotData(hSFC(arraySide),freqValsMT,dataSFC,colorNamesAttCue(attCuePos,:));
                plotData(hAmpCorr(arraySide),freqVals,dataAmpCorr,colorNamesAttCue(attCuePos,:));
                plotData(hffPPC(arraySide),freqValsPPC,dataffPPC,colorNamesAttCue(attCuePos,:));
                plotData(hsfPPC(arraySide),freqValsPPC,datasfPPC,colorNamesAttCue(attCuePos,:));                
                
            end
            
            makeBarPlot(hBarRsc(arraySide),rSCData{arraySide},colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarFFCAlpha(arraySide),mean(ffcData{arraySide}(:,:,alphaPos),3),colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarFFCSSVEP(arraySide),ffcData{arraySide}(:,:,ssvepPos),colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarSFCAlpha(arraySide),mean(sfcData{arraySide}(:,:,alphaPos),3),colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarSFCSSVEP(arraySide),sfcData{arraySide}(:,:,ssvepPos),colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarAmpCorrAlpha(arraySide),mean(ampCorrData{arraySide}(:,:,alphaPos),3),colorNamesAttCue,showAbsoluteValsFlag);
            makeBarPlot(hBarAmpCorrSSVEP(arraySide),ampCorrData{arraySide}(:,:,ssvepPos),colorNamesAttCue,showAbsoluteValsFlag);
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
        
        yLims = getYLims(hFFC);
        axis(hFFC(1),[str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String')) yLims]);
        if showAbsoluteValsFlag;    set(hFFC(1),'YLim',[0 1]);          end
        
        yLims = getYLims(hSFC);
        axis(hSFC(1),[str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String')) yLims]);
        if showAbsoluteValsFlag;    set(hSFC(1),'YLim',[0 1]);          end
        
        yLims = getYLims(hAmpCorr);
        axis(hAmpCorr(1),[str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String')) yLims]);
        if showAbsoluteValsFlag;    set(hAmpCorr(1),'YLim',[0 1]);      end
        
        yLims = getYLims(hffPPC);
        axis(hffPPC(1),[str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String')) yLims]);
        if showAbsoluteValsFlag;    set(hffPPC(1),'YLim',[0 1]);          end
        
        yLims = getYLims(hsfPPC);
        axis(hsfPPC(1),[str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String')) yLims]);
        if showAbsoluteValsFlag;    set(hsfPPC(1),'YLim',[0 1]);          end        

        yLims = getYLims(hBarFR(1:2)); axis(hBarFR(1),[0 6 yLims]);
        yLims = getYLims(hBarAlpha(1:2)); axis(hBarAlpha(1),[0 6 yLims]);
        yLims = getYLims(hBarSSVEP(1:2)); axis(hBarSSVEP(1),[0 6 yLims]);
        ylim(hBehavior(1),[0 1]);
        
        if showAbsoluteValsFlag
            yLims = [0 1];
        else
            yLims = [-0.05 0.05];
        end
        axis(hBarRsc(1),[0 6 yLims]);
        axis(hBarFFCAlpha(1),[0 6 yLims]);
        axis(hBarFFCSSVEP(1),[0 6 yLims]);
        axis(hBarSFCAlpha(1),[0 6 yLims]);
        axis(hBarSFCSSVEP(1),[0 6 yLims]);
        axis(hBarAmpCorrAlpha(1),[0 6 yLims]);
        axis(hBarAmpCorrSSVEP(1),[0 6 yLims]);
        
        % Labels
        xlabel(hPSTH(3),'Time (s)'); title(hPSTH(1),'Firing Rate (spikes/s)');
        xlabel(hERP(3),'Time (s)'); title(hERP(1),'ERP (\muV)');
        xlabel(hFFT(3),'Frequency (Hz)'); title(hFFT(1),'Log FFT');
        
        title(hBarFR(1),'Firing Rate');
        title(hBarAlpha(1),'Alpha');
        title(hBarSSVEP(1),'SSVEP');
        
        xlabel(hFFC(4),'Frequency (Hz)'); title(hFFC(1),'FFC');
        xlabel(hSFC(4),'Frequency (Hz)'); title(hSFC(1),'SFC');
        
        xlabel(hAmpCorr(4),'Frequency (Hz)'); title(hAmpCorr(1),'Amp Corr');
        xlabel(hffPPC(4),'Frequency (Hz)'); title(hffPPC(1),'ffPPC');
        xlabel(hsfPPC(4),'Frequency (Hz)'); title(hsfPPC(1),'sfPPC');
        
        title(hBarRsc(1),'r_s_c');
        title(hBarFFCAlpha(1),'FFC_A_l_p_h_a');
        title(hBarFFCSSVEP(1),'FFC_S_S_V_E_P');
        title(hBarSFCAlpha(1),'SFC_A_l_p_h_a');
        title(hBarSFCSSVEP(1),'SFC_S_S_V_E_P');
        title(hBarAmpCorrAlpha(1),'amp_A_l_p_h_a');
        title(hBarAmpCorrSSVEP(1),'amp_S_S_V_E_P');

        title(hBehavior(1),'Behavior (Fraction Correct)');
        
        for arraySide=1:3
            ylabel(hPSTH(arraySide),['N=' num2str(size(psthData{arraySide},2))],'color',colorNamesSides(arraySide));
        end
        
        for arraySide = 1:4
            ylabel(hFFC(arraySide),['N=' num2str(size(ffcData{arraySide},2))],'color',colorNamesSides(arraySide));
        end
        
        legend(hBehavior(1),'0V','1V','0I','1I','N','location','southeast');
    end
    function rescaleXY_Callback(~,~)
        tRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        fRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        
        axis(hPSTH(1),[tRange str2double(get(hFRRangeMin,'String')) str2double(get(hFRRangeMax,'String'))]);
        axis(hERP(1),[tRange str2double(get(hERPRangeMin,'String')) str2double(get(hERPRangeMax,'String'))]);
        axis(hFFT(1),[fRange str2double(get(hFFTYMin,'String')) str2double(get(hFFTYMax,'String'))]);
        
        xlim(hFFC(1),fRange);
        xlim(hSFC(1),fRange);
        xlim(hAmpCorr(1),fRange);
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
        
        holdOnGivenPlotHandle(hFFC,holdOnState);
        holdOnGivenPlotHandle(hSFC,holdOnState);
        holdOnGivenPlotHandle(hffPPC,holdOnState);
        holdOnGivenPlotHandle(hsfPPC,holdOnState);
        holdOnGivenPlotHandle(hAmpCorr,holdOnState);
        holdOnGivenPlotHandle(hBarRsc,holdOnState);
        holdOnGivenPlotHandle(hBarFFCAlpha,holdOnState);
        holdOnGivenPlotHandle(hBarFFCSSVEP,holdOnState);
        holdOnGivenPlotHandle(hBarSFCAlpha,holdOnState);
        holdOnGivenPlotHandle(hBarSFCSSVEP,holdOnState);
        holdOnGivenPlotHandle(hBarAmpCorrAlpha,holdOnState);
        holdOnGivenPlotHandle(hBarAmpCorrSSVEP,holdOnState);
        
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
        
        claGivenPlotHandle(hFFC);
        claGivenPlotHandle(hSFC);
        claGivenPlotHandle(hffPPC);
        claGivenPlotHandle(hsfPPC);
        claGivenPlotHandle(hAmpCorr);
        claGivenPlotHandle(hBarRsc);
        claGivenPlotHandle(hBarFFCAlpha);
        claGivenPlotHandle(hBarFFCSSVEP);
        claGivenPlotHandle(hBarSFCAlpha);
        claGivenPlotHandle(hBarSFCSSVEP);
        claGivenPlotHandle(hBarAmpCorrAlpha);
        claGivenPlotHandle(hBarAmpCorrSSVEP);
        
        claGivenPlotHandle(hElectrodes);
        
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
function [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ampCorrData,rSCData,ffcData,ffPhiData,sfcData,sfPhiData,freqValsMT,ffppcData,sfppcData,freqValsPPC,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getData(folderSourceString,fileNameStringTMP,neuronType,populationType,tStr,oStr,tpStr,tapers)

numDatasets = length(fileNameStringTMP);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ampCorrData,rSCData,ffcData,ffPhiData,sfcData,sfPhiData,freqValsMT,ffppcData,sfppcData,freqValsPPC,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getDataSingleSession(folderSourceString,fileNameStringTMP{1},neuronType,populationType,tStr,oStr,tpStr,tapers); % First session

if length(fileNameStringTMP)>1
    for i=2:numDatasets
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [psthDataTMP,xsFRTMP,firingRatesTMP,erpDataTMP,timeValsTMP,fftDataTMP,freqValsTMP,alphaDataTMP,ssvepDataTMP,ampCorrDataTMP,rSCDataTMP,ffcDataTMP,ffPhiDataTMP,sfcDataTMP,sfPhiDataTMP,freqValsMTTMP,ffppcDataTMP,sfppcDataTMP,freqValsPPCTMP,electrodeArrayTMP,perCorrectTMP,uniqueOrientationChangeDegTMP] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},neuronType,populationType,tStr,oStr,tpStr,tapers);
        
        perCorrect = perCorrect + perCorrectTMP;
        uniqueOrientationChangeDeg = uniqueOrientationChangeDeg + uniqueOrientationChangeDegTMP;
        
        for k=1:2
            electrodeArray{k} = cat(1,electrodeArray{k},electrodeArrayTMP{k});
        end
        
        for k=1:3 % for each array side and 3- both arrays combined
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
            else
                error('freqVals do not match');
            end
        end
        
        for k=1:4 % 1- Right Array elec pairs 2-Left Array elec pairs 3-Intra-hemispheric pairwise combined 4-inter-hemispheric pairwise rSC
            if isequal(freqValsMT,freqValsMTTMP)
                ffcData{k} = cat(2,ffcData{k},ffcDataTMP{k});
                ffPhiData{k} = cat(2,ffPhiData{k},ffPhiDataTMP{k});
                sfcData{k} = cat(2,sfcData{k},sfcDataTMP{k});
                sfPhiData{k} = cat(2,sfPhiData{k},sfPhiDataTMP{k});
            else
                error('freqValsMT do not match')
            end
            
            if isequal(freqVals,freqValsTMP)
                ampCorrData{k} = cat(2,ampCorrData{k},ampCorrDataTMP{k});
            else
                error('freqVals do not match');
            end
            
            if isequal(freqValsPPC,freqValsPPCTMP)
                ffppcData{k} = cat(2,ffppcData{k},ffppcDataTMP{k});
                sfppcData{k} = cat(2,sfppcData{k},sfppcDataTMP{k});
            else
                error('freqValsPPC do not match')
            end
            
            rSCData{k} = cat(2,rSCData{k},rSCDataTMP{k});
        end
    end
    perCorrect = perCorrect/numDatasets;
    uniqueOrientationChangeDeg = uniqueOrientationChangeDeg/numDatasets;
end
end
function [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,ampCorrData,rSCData,ffcData,ffPhiData,sfcData,sfPhiData,freqValsMT,ffppcData,sfppcData,freqValsPPC,electrodeArray,perCorrect,uniqueOrientationChangeDeg] = getDataSingleSession(folderSourceString,fileNameString,neuronType,populationType,tStr,oStr,tpStr,tapers)

folderSave = fullfile(folderSourceString,'Data','savedData');
makeDirectory(folderSave);

fileToSave = fullfile(folderSave,[fileNameString neuronType populationType tStr oStr tpStr '.mat']);
coherencyFileToSave = fullfile(folderSave,[fileNameString neuronType populationType tStr oStr tpStr '_coherenceData_tapers_',num2str(tapers(1)) '.mat']);

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
    [electrodepairsWithinHemisphere,electrodePairsAcrossHemispheres] = getGoodElectrodePairs(electrodeArray);
    
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
        elseif strcmp(tpStr(1:11),'TargetOnset')
            if strcmp(tpStr,'TargetOnset_250ms')
                timeRange = [-0.25 0];
            elseif strcmp(tpStr,'TargetOnset_500ms')
                timeRange = [-0.5 0];
            end  
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
    
    for i=1:2 % Each array side
        eList = electrodeArray{i};
        ePairList = electrodepairsWithinHemisphere{i};
        clear psthDataTMP firingRatesTMP erpDataTMP fftDataTMP alphaDataTMP ssvepDataTMP ampCorrTMP rSCTMP
        clear ffcTMP ffPhiTMP sfcTMP sfPhiTMP ffppcTMP sfppcTMP
        
        for j=1:numConditions
            disp(['numArray:' num2str(i) ', numCondition:' num2str(j)]);
            
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
            
            % Amplitude Correlation
            ampCorrTMP(j,:,:) = getAmplitudeCorrelation(lfpData{j}.segmentedLFPData(:,:,pos),ePairList);
            
            % Spike-count Correlation (rSC)
            rSCTMP(j,:) = getSpikeCountCorrelation(spikeData{j}.segmentedSpikeData,ePairList,timeRange);
            
            % Coherency and Spike-LFP Phase Analysis Within Hemisphere(FFC,SFC,sfPhi)
            [ffcTMP(j,:,:),ffPhiTMP(j,:,:),sfcTMP(j,:,:),sfPhiTMP(j,:,:),freqValsMT] = getCoherencyMeasures(lfpData{j}.segmentedLFPData(:,:,pos),spikeData{j}.segmentedSpikeData,ePairList,tapers,lfpData{j}.timeVals(pos),timeRange);
            
            % Pairwise Phase consistency(FF-PPC & SF-PPC analysis within Hemisphere) 
            [ffppcTMP(j,:,:),sfppcTMP(j,:,:),freqValsPPC] = getPairWisePhaseConsistencyMeasures(lfpData{j}.segmentedLFPData(:,:,pos),spikeData{j}.segmentedSpikeData,ePairList,tapers,lfpData{j}.timeVals(pos),timeRange);

        end
        
        psthData{i} = psthDataTMP;
        firingRates{i} = firingRatesTMP;
        erpData{i} = erpDataTMP;
        fftData{i} = fftDataTMP;
        alphaData{i} = alphaDataTMP;
        ssvepData{i} = ssvepDataTMP;
        ampCorrData{i} = ampCorrTMP;
        rSCData{i} = rSCTMP;
        
        ffcData{i} = ffcTMP;
        ffPhiData{i} = ffPhiTMP;
        sfcData{i} = sfcTMP;
        sfPhiData{i} = sfPhiTMP;
        
        ffppcData{i} = ffppcTMP;
        sfppcData{i} = sfppcTMP;
    end
    
    % Combine Data across both arrays/hemispheres
    combinedDataPos = 3;
    psthData{combinedDataPos} = combineDataAcrossBothArrays(psthData);
    firingRates{combinedDataPos} = combineDataAcrossBothArrays(firingRates);
    erpData{combinedDataPos} = combineDataAcrossBothArrays(erpData);
    fftData{combinedDataPos} = combineDataAcrossBothArrays(fftData);
    alphaData{combinedDataPos} = combineDataAcrossBothArrays(alphaData);
    ssvepData{combinedDataPos} = combineDataAcrossBothArrays(ssvepData);
    ampCorrData{combinedDataPos} = combineDataAcrossBothArrays(ampCorrData);
    rSCData{combinedDataPos} = combineDataAcrossBothArrays(rSCData);
    
    ffcData{combinedDataPos} = combineDataAcrossBothArrays(ffcData);
    ffPhiData{combinedDataPos} = combineDataAcrossBothArrays(ffPhiData);
    sfcData{combinedDataPos} = combineDataAcrossBothArrays(sfcData);
    sfPhiData{combinedDataPos} = combineDataAcrossBothArrays(sfPhiData);
    
    ffppcData{combinedDataPos} = combineDataAcrossBothArrays(ffppcData);
    sfppcData{combinedDataPos} = combineDataAcrossBothArrays(sfppcData);
    
    % AmpCorr, rSC & Coherency analysis across hemispheres
    disp('Working on amplitude correlation Data, rSC Data and Coherency Data for electrode pairs across hemispheres')
    for j = 1:numConditions
        % Amplitude Correlation analysis across Hemispheres
        ampCorrData_AH(j,:,:) = getAmplitudeCorrelation(lfpData{j}.segmentedLFPData(:,:,pos),electrodePairsAcrossHemispheres);
        % rSC Analysis Across Hemispheres
        rSCData_AH(j,:) = getSpikeCountCorrelation(spikeData{j}.segmentedSpikeData,electrodePairsAcrossHemispheres,timeRange);
        % Coherency and Spike-LFP Phase Analysis Across Hemispheres(FFC,SFC,sfPhi)
        [ffcData_AH(j,:,:),ffPhiData_AH(j,:,:),sfcData_AH(j,:,:),sfPhiData_AH(j,:,:),~] = getCoherencyMeasures(lfpData{j}.segmentedLFPData(:,:,pos),spikeData{j}.segmentedSpikeData,electrodePairsAcrossHemispheres,tapers,lfpData{j}.timeVals(pos),timeRange); %AH -across Hemisheres
        % PairWise Phase Consistency ffPPC and sfPPC across Hemispheres
        [ffppcData_AH(j,:,:),sfppcData_AH(j,:,:)] = getPairWisePhaseConsistencyMeasures(lfpData{j}.segmentedLFPData(:,:,pos),spikeData{j}.segmentedSpikeData,electrodePairsAcrossHemispheres,tapers,lfpData{j}.timeVals(pos),timeRange);
        
    end
    
    % Combining Coherency Data within Hemispheres and across hemispheres together
    InterHemisphericCoherencyPos = 4;
    ffcData{InterHemisphericCoherencyPos} = ffcData_AH;
    ffPhiData{InterHemisphericCoherencyPos} = ffPhiData_AH;
    sfcData{InterHemisphericCoherencyPos} = sfcData_AH;
    sfPhiData{InterHemisphericCoherencyPos} = sfPhiData_AH;
    rSCData{InterHemisphericCoherencyPos} = rSCData_AH;
    ampCorrData{InterHemisphericCoherencyPos} = ampCorrData_AH;
    ffppcData{InterHemisphericCoherencyPos} = ffppcData_AH;
    sfppcData{InterHemisphericCoherencyPos} = sfppcData_AH;
    
    % Save data
    save(fileToSave,'psthData','xsFR','firingRates','erpData','timeVals','fftData','freqVals','alphaData','ssvepData','ampCorrData','rSCData','electrodeArray','perCorrect','uniqueOrientationChangeDeg');
    save(coherencyFileToSave,'ffcData','ffPhiData','sfcData','sfPhiData','freqValsMT','ffppcData','sfppcData','freqValsPPC');
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
function [ffc,ffPhi,sfc,sfPhi,freqValsMT] = getCoherencyMeasures(lfpData,spikeData,electrodePair,tapers,timeVals,timeRange)

% Set up MT
Fs              = round(1/(timeVals(2)-timeVals(1)));
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 1;

disp('Working on FFC Data')
% Field-Field coherence
for i=1:size(electrodePair,1)
    clear lfp1 lfp2
%    disp(['FFC ElectrodePair:',num2str(i)]);
    lfp1 = squeeze(lfpData(electrodePair(i,1),:,:));
    lfp2 = squeeze(lfpData(electrodePair(i,2),:,:));
    
    [ffc(i,:),ffPhi(i,:),~,~,~,freqFFC]=coherencyc(lfp1',lfp2',params); %#ok<*AGROW>
end

% if params.trialave = 0;
% [ppc_ffc,PLF_ffc,without_phaseCovariation_ffc,traditional_ffc]= computeCoherencyFromSpectrum(S1_ffc,S2_ffc,S12_ffc,params);
% ffc = traditional_ffc;
% end

disp('Working on SFC Data')
% Spike-Field coherence
for i=1:size(electrodePair,1)
    clear lfp spk
%    disp(['SFC ElectrodePair:',num2str(i)]);
    lfp = squeeze(lfpData(electrodePair(i,1),:,:));
    spk = convertSpikeTimes2Bins(spikeData(electrodePair(i,2),:,:),timeRange,1000/Fs);
    [sfc(i,:),sfPhi(i,:),~,~,~,freqSFC]=coherencycpb(lfp',spk,params); %#ok<*AGROW>
end

% Sanity Check
if isequal(freqFFC,freqSFC)
    freqValsMT = freqFFC;
else
    error('freqVals from FFC & SFC do not match!')
end

end
function [ffPPC,sfPPC,freqValsPPC] = getPairWisePhaseConsistencyMeasures(lfpData,spikeData,electrodePair,tapers,timeVals,timeRange)

% PPC- adapted from fieldtrip/connectivity/ft_connectivity_ppc.m computes 
% pairwise phase consistency  from a data-matrix containing a cross-spectral 
% density. This implements the method described in Vinck M, van Wingerden M, 
% Womelsdorf T, Fries P, Pennartz CM.
% The pairwise phase consistency: a bias-free measure of rhythmic neuronal
% synchronization. Vinck et al. Neuroimage. 2010

% Set up MT
Fs              = round(1/(timeVals(2)-timeVals(1)));
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 0;

disp('Working on FF-PPC')
% Field-Field Pairwise Phase Consistency
for i=1:size(electrodePair,1)
    clear lfp1 lfp2 input 
%    disp(['ffPPC ElectrodePair:',num2str(i)]);
    lfp1 = squeeze(lfpData(electrodePair(i,1),:,:));
    lfp2 = squeeze(lfpData(electrodePair(i,2),:,:));
    
    [~,~,S12_ffc,~,~,freqValsFFPPC]=coherencyc(lfp1',lfp2',params); %#ok<*AGROW>
    %ffPPC
    input = (S12_ffc./abs(S12_ffc))'; % normalize the cross-spectrum
    siz = size(input);
    n = siz(1);
    if n>1
        outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
        ffPPC(i,:)  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
    else
        error('computation of PPC requires >1 trial, please feed all trial dataset into computeCoherencyFromSpectrum program')
    end
end

disp('Working on SF-PPC')
% Spike-Field Pairwise Phase Consistency
for i=1:size(electrodePair,1)
    clear lfp spk input 
%    disp(['sfPPC ElectrodePair:',num2str(i)]);
    lfp = squeeze(lfpData(electrodePair(i,1),:,:));
    spk = convertSpikeTimes2Bins(spikeData(electrodePair(i,2),:,:),timeRange,1000/Fs);
    [~,~,S12_sfc,~,~,freqValsSFPPC]=coherencycpb(lfp',spk,params); %#ok<*AGROW>
    %sfPPC
    input = (S12_sfc./abs(S12_sfc))'; % normalize the cross-spectrum
    siz = size(input);
    n = siz(1);
    if n>1
        outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
        sfPPC(i,:)  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
    else
        error('computation of PPC requires >1 trial, please feed all trial dataset into computeCoherencyFromSpectrum program')
    end
    
end

% Sanity Check
if isequal(freqValsFFPPC,freqValsSFPPC)
    freqValsPPC = freqValsFFPPC;
else
    error('freqVals from FFC & SFC do not match!')
end

end
function ampCorr = getAmplitudeCorrelation(lfpData,electrodePair)
disp('Working on amplitude Correlation Data')
fftDataTrialWise = abs(fft(lfpData,[],3));
for k = 1:size(electrodePair,1)
    clear fft1 fft2
%    disp(['AmpCorrData,electrodePair:',num2str(k)])
%    disp(['Amplitude Correlation ElectrodePair:',num2str(k)])
    fft1 = squeeze(fftDataTrialWise(electrodePair(k,1),:,:));
    fft2 = squeeze(fftDataTrialWise(electrodePair(k,2),:,:));
    ampCorr(k,:) = diag(corr(fft1,fft2))';
end
end
function rSC = getSpikeCountCorrelation(spikeData,electrodePair,tRangeS) % Computes for all Electrode pairs
disp('Working on rSc Data')
for k=1:size(electrodePair,1)
    clear h1 h2
    h1 = getSpikeCounts(spikeData(electrodePair(k,1),:),tRangeS);
    h2 = getSpikeCounts(spikeData(electrodePair(k,2),:),tRangeS);
    rSC(k) = (mean(h1.*h2) - mean(h1)*mean(h2))/(std(h1)*std(h2));
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
function makeBarPlot(h,data,colorNames,showAbsoluteValsFlag)

N = size(data,2);
if ~showAbsoluteValsFlag
    data = data - repmat(data(5,:),5,1);
end

mData = mean(data,2);
semData = std(data,[],2)/sqrt(N);

for i=1:size(data,1)
    if i==3 || i==4
        continue
    end
    plot(h,i,mData(i),'color',colorNames(i,:),'marker','o');
    hold(h,'on');
    errorbar(h,i,mData(i),semData(i),'color',colorNames(i,:));
end
set(h,'XTick',[],'XTicklabel',[]);
end
function plotData(hPlot,xs,data,colorName)

if isequal(colorName,[0 0 1])
    colorName2 = [0 1 1];
elseif isequal(colorName,[1 0 0])
    colorName2 = [1 0 1];
else
    colorName2 = colorName;
end

mData = mean(data,1);
sData = std(data,[],1)/sqrt(size(data,1));
xsLong = [xs fliplr(xs)];
ysLong = [mData+sData fliplr(mData-sData)];

patch(xsLong,ysLong,colorName2,'parent',hPlot);
hold(hPlot,'on');
plot(hPlot,xs,mData,'color',colorName,'linewidth',2); 

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