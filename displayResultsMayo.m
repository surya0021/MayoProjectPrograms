function displayResultsMayo(folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.8; backgroundColor = 'w'; % Panels

% Show electrodes
electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
hElectrodes = showElectrodeLocationsMayo(electrodeGridPos,[],'r',[],0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hParameterPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.25 panelStartHeight 0.25 panelHeight]);
paramsHeight=1/7;

% FileNameString
[fileNameStringAll,fileNameStringListArray] = getFileNameStringList;
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

% CueType
cueTypeString = [{'Valid'} {'Invalid'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-5*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Cue Type','FontSize',fontSizeSmall);
hCueType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [0.5 1-5*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',cueTypeString,'FontSize',fontSizeSmall);

% Analysis Interval
analysisIntervalString = [{'Baseline'} {'StimOnset'} {'TargetOnset'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-6*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);
hTimePeriodType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-6*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',analysisIntervalString ,'FontSize',fontSizeSmall);

% Orientation Change
orientationChangeString=[{'Selected'} {'All'}];
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 1-7*paramsHeight 0.5 paramsHeight], ...
    'Style','text','String','Orientation Change','FontSize',fontSizeSmall);
hOriChangeType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 1-7*paramsHeight 0.5 paramsHeight], ...
    'Style','popup','String',orientationChangeString ,'FontSize',fontSizeSmall);

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
hBehavior = getPlotHandles(3,1,[0.05 0.05 0.2 0.7],0,0.05,0); linkaxes(hBehavior);

hPSTH = getPlotHandles(3,1,[0.3 0.05 0.1 0.7],0,0.05,0); linkaxes(hPSTH);
hERP = getPlotHandles(3,1,[0.45 0.05 0.1 0.7],0,0.05,0); linkaxes(hERP);
hFFT = getPlotHandles(3,1,[0.6 0.05 0.1 0.7],0,0.05,0); linkaxes(hFFT);

hBarFR = getPlotHandles(3,1,[0.75 0.05 0.05 0.7],0,0.05,0); linkaxes(hBarFR);
hBarAlpha = getPlotHandles(3,1,[0.82 0.05 0.05 0.7],0,0.05,0); linkaxes(hBarAlpha);
hBarSSVEP = getPlotHandles(3,1,[0.9 0.05 0.05 0.7],0,0.05,0); linkaxes(hBarSSVEP);

colorNamesSides = 'cm';

% Plotting functions
    function plotData_Callback(~,~)
        
        s=get(hSession,'val'); fileNameStringTMP = fileNameStringListArray{s};
        n=get(hNeuronType,'val'); neuronType = neuronTypeString{n};
        p=get(hPopulationType,'val'); populationType = populationString{p};
        t=get(hTrialOutcome,'val'); tStr = trialOutcomeString{t}(1);
        c=get(hCueType,'val'); cStr = cueTypeString{c}(1);
        tp=get(hTimePeriodType,'val'); tpStr = analysisIntervalString{tp};
        o=get(hOriChangeType,'val'); oStr=orientationChangeString{o};
        colorNamesAttPos = colorNames{get(hChooseColor,'val')};
        holdState=get(hHoldOn,'Value');
        if strcmp(tpStr,'TargetOnset')
            tRange = [-0.5 0.1];
        else
            tRange = [-0.25 0.5];
        end
        
        % Get data
        [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,electrodeArray,perCorrect] = getData(folderSourceString,fileNameStringTMP,neuronType,populationType,tStr,cStr,tpStr,oStr);
        if holdState
            x=0.4;
        else
            x=0.05;
        end
        for arraySide=1:2
            showElectrodeLocationsMayo([],electrodeArray{arraySide},colorNamesSides(arraySide),hElectrodes,1,0); % Show electrodes used for analysis
            
            for attPos=1:3
                plot(hPSTH(arraySide),xsFR,squeeze(mean(psthData{arraySide}(attPos,:,:),2)),'color',colorNamesAttPos(attPos));
                hold(hPSTH(arraySide),'on');
                
                plot(hERP(arraySide),timeVals,squeeze(mean(erpData{arraySide}(attPos,:,:),2)),'color',colorNamesAttPos(attPos)); 
                hold(hERP(arraySide),'on');
                
                plot(hFFT(arraySide),freqVals,squeeze(mean(fftData{arraySide}(attPos,:,:),2)),'color',colorNamesAttPos(attPos)); 
                hold(hFFT(arraySide),'on');
            
                if arraySide==1
                      if attPos==3;
                        text(x,0.35-0.1*attPos,[tStr 'N'],'Unit','normalized','color',colorNamesAttPos(attPos),'parent',hPSTH(1));
                        text(x,0.35-0.1*attPos,'Neutral','Unit','normalized','color',colorNamesAttPos(attPos),'parent',hPSTH(3));
                        plot(hPSTH(3),xsFR,squeeze(mean(cat(2,psthData{1}(attPos,:,:),psthData{2}(attPos,:,:)),2)),'color',colorNamesAttPos(attPos));
                        hold(hPSTH(3),'on');
                        
                        plot(hERP(3),timeVals,squeeze(mean(cat(2,erpData{1}(attPos,:,:),erpData{2}(attPos,:,:)),2)),'color',colorNamesAttPos(attPos));
                        hold(hERP(3),'on');
                        
                        plot(hFFT(3),freqVals,squeeze(mean(cat(2,fftData{1}(attPos,:,:),fftData{2}(attPos,:,:)),2)),'color',colorNamesAttPos(attPos));
                        hold(hFFT(3),'on');
                      else
                        str=[{'In'} {'Out'}];
                        text(x,0.35-0.1*attPos,[tStr num2str(attPos-1) cStr],'Unit','normalized','color',colorNamesAttPos(attPos),'parent',hPSTH(1)); 
                       
                        plot(hPSTH(3),xsFR,squeeze(mean(cat(2,psthData{1}(attPos,:,:),psthData{2}(3-attPos,:,:)),2)),'color',colorNamesAttPos(attPos));
                        text(x,0.35-0.1*attPos,['Att ' str{attPos} ' ' cStr],'Unit','normalized','color',colorNamesAttPos(attPos),'parent',hPSTH(3)); 
                        hold(hPSTH(3),'on');
                        
                        plot(hERP(3),timeVals,squeeze(mean(cat(2,erpData{1}(attPos,:,:),erpData{2}(3-attPos,:,:)),2)),'color',colorNamesAttPos(attPos));
                        hold(hERP(3),'on');
                        
                        plot(hFFT(3),freqVals,squeeze(mean(cat(2,fftData{1}(attPos,:,:),fftData{2}(3-attPos,:,:)),2)),'color',colorNamesAttPos(attPos));
                        hold(hFFT(3),'on');
                      end
                end
            end
            makeBarPlot(hBarFR(arraySide),firingRates{arraySide},colorNamesAttPos);
            makeBarPlot(hBarAlpha(arraySide),alphaData{arraySide},colorNamesAttPos);
            makeBarPlot(hBarSSVEP(arraySide),ssvepData{arraySide},colorNamesAttPos);
            
            if isempty(get(hBehavior(arraySide),'Children'))
                makeBehaviorPlot(hBehavior(arraySide),perCorrect{arraySide},colorNamesAttPos);
            end
        end
       
        for i=1:3
            if i==3
                firingRatesBothArrays(i,:)=cat(2,firingRates{1}(i,:),firingRates{2}(i,:));
                alphaDataBothArrays(i,:)=cat(2,alphaData{1}(i,:),alphaData{2}(i,:));
                ssvepDataBothArrays(i,:)=cat(2,ssvepData{1}(i,:),ssvepData{2}(i,:));
            else
                firingRatesBothArrays(i,:)=cat(2,firingRates{1}(i,:),firingRates{2}(3-i,:));
                alphaDataBothArrays(i,:)=cat(2,alphaData{1}(i,:),alphaData{2}(3-i,:));
                ssvepDataBothArrays(i,:)=cat(2,ssvepData{1}(i,:),ssvepData{2}(3-i,:));
            end
        end
        
         makeBarPlot(hBarFR(3),firingRatesBothArrays,colorNamesAttPos);
         makeBarPlot(hBarAlpha(3),alphaDataBothArrays,colorNamesAttPos);
         makeBarPlot(hBarSSVEP(3),ssvepDataBothArrays,colorNamesAttPos);
         
         if isempty(get(hBehavior(3),'Children'))
            makeBehaviorPlot(hBehavior(3),perCorrect{3},colorNamesAttPos);
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
        
        yLims = getYLims(hBarFR(1:2)); axis(hBarFR(1),[0 4 0 yLims(2)]);
        yLims = getYLims(hBarAlpha(1:2)); axis(hBarAlpha(1),[0 4 yLims]);
        yLims = getYLims(hBarSSVEP(1:2)); axis(hBarSSVEP(1),[0 4 yLims]);
        
        % Labels
        xlabel(hPSTH(3),'Time (s)'); title(hPSTH(1),'Firing Rate (spikes/s)');
        xlabel(hERP(3),'Time (s)'); title(hERP(1),'ERP (\muV)');
        xlabel(hFFT(3),'Frequency (Hz)'); title(hFFT(1),'Log FFT');
        xlabel(hBehavior(3),'Orientation Change Level'); title(hBehavior(1),'Detection performance (%Correct)');
        
        title(hBarFR(1),'Firing Rate');
        title(hBarAlpha(1),'Alpha');
        title(hBarSSVEP(1),'SSVEP');
        
        for arraySide=1:2
            ylabel(hPSTH(arraySide),['N=' num2str(length(electrodeArray{arraySide}))],'color',colorNamesSides(arraySide));
        end
        
        ylabel(hBehavior(1),'Right Array','FontWeight','bold');
        ylabel(hBehavior(2),'Left Array','FontWeight','bold');
        ylabel(hBehavior(3),'Both Arrays','FontWeight','bold');
        ylabel(hPSTH(3),['N=' num2str(length(electrodeArray{1})+length(electrodeArray{2}))],'color','b');
        %Legends
        legend(hBehavior(1),'Cued','Uncued','Neutral','Location','southeast');
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
        showElectrodeLocationsMayo([],1:96,'w',hElectrodes,1,0);
        legend(hBehavior(1),'off');
        for n=1:3
        ylabel(hPSTH(n),'');
        end
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

function [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,electrodeArray,perCorrect] = getData(folderSourceString,fileNameStringTMP,neuronType,populationType,tStr,cStr,tpStr,oStr)

disp(['Working on dataset 1 of ' num2str(length(fileNameStringTMP))]);
[psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,electrodeArray,perCorrect] = getDataSingleSession(folderSourceString,fileNameStringTMP{1},neuronType,populationType,tStr,cStr,tpStr,oStr); % First session

if length(fileNameStringTMP)>1
    for i=2:length(fileNameStringTMP)
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [psthDataTMP,xsFRTMP,firingRatesTMP,erpDataTMP,timeValsTMP,fftDataTMP,freqValsTMP,alphaDataTMP,ssvepDataTMP,electrodeArrayTMP,perCorrectTMP] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},neuronType,populationType,tStr,cStr,tpStr,oStr);
        
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
            else
                error('freqVals do not match');
            end
            
            electrodeArray{k} = cat(1,electrodeArray{k},electrodeArrayTMP{k});
            perCorrect{k}=cat(3,perCorrect{k},perCorrectTMP{k});
        end
            perCorrect{3}=cat(3,perCorrect{3},perCorrectTMP{3});
    end
end
end
function [psthData,xsFR,firingRates,erpData,timeVals,fftData,freqVals,alphaData,ssvepData,electrodeArray,perCorrect] = getDataSingleSession(folderSourceString,fileNameString,neuronType,populationType,tStr,cStr,tpStr,oStr)

binWidthMS = 10;
alphaRangeHz = [8 12]; ssvepFreqHz = 20;

folderName = fullfile(folderSourceString,'Data','segmentedData',fileNameString);

electrodeArray = getGoodElectrodes(folderSourceString,fileNameString,neuronType,populationType);

% Get Spike and LFP data
for attLoc=0:1
    
    if strcmp(tpStr,'Baseline')
        timeRange = [-0.25 0];
        lfpData{attLoc+1}=load(fullfile(folderName,[fileNameString tStr num2str(attLoc) cStr '_StimOnset_LFP']));
        spikeData{attLoc+1}=load(fullfile(folderName,[fileNameString tStr num2str(attLoc) cStr '_StimOnset_Spikes']));
    elseif strcmp(tpStr,'StimOnset')
        timeRange = [0.25 0.5];
        lfpData{attLoc+1}=load(fullfile(folderName,[fileNameString tStr num2str(attLoc) cStr '_StimOnset_LFP']));
        spikeData{attLoc+1}=load(fullfile(folderName,[fileNameString tStr num2str(attLoc) cStr '_StimOnset_Spikes']));
    elseif strcmp(tpStr,'TargetOnset')
        timeRange = [-0.5 0];
        lfpData{attLoc+1}=load(fullfile(folderName,[fileNameString tStr num2str(attLoc) cStr '_TargetOnset_LFP']));
        spikeData{attLoc+1}=load(fullfile(folderName,[fileNameString tStr num2str(attLoc) cStr '_TargetOnset_Spikes']));
    end
    if strcmp(oStr,'Selected')
        goodPos{attLoc+1}=getGoodTrials(fileNameString,[tStr num2str(attLoc) cStr]);
    else
        goodPos{attLoc+1}=(1:size(lfpData{attLoc+1}.segmentedLFPData,2));
    end
end

if strcmp(tpStr,'TargetOnset')
    lfpData{3}=load(fullfile(folderName,[fileNameString tStr 'N_TargetOnset_LFP']));
    spikeData{3}=load(fullfile(folderName,[fileNameString tStr 'N_TargetOnset_Spikes']));
else
    lfpData{3}=load(fullfile(folderName,[fileNameString tStr 'N_StimOnset_LFP']));
    spikeData{3}=load(fullfile(folderName,[fileNameString tStr 'N_StimOnset_Spikes']));
end
if strcmp(oStr,'Selected')
    goodPos{3}=getGoodTrials(fileNameString,[tStr 'N']);
else
    goodPos{3}=(1:size(lfpData{3}.segmentedLFPData,2));
end
for i=1:3
    disp([fileNameString tStr num2str(i-1) cStr ', Stimulus repeats for Loc' num2str(i-1) ': ' num2str(size(goodPos{i},2))]);
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
    
    clear psthDataTMP firingRatesTMP erpDataTMP fftDataTMP alphaDataTMP ssvepDataTMP
    for attPos=1:3
        % Firing Rates
        for k=1:length(eList)
            [psthDataTMP(attPos,k,:),xsFR] = getPSTH(spikeData{attPos}.segmentedSpikeData(eList(k),goodPos{attPos}),binWidthMS,[timeVals(1) timeVals(end)]);
            firingRatesTMP(attPos,k) = mean(getSpikeCounts(spikeData{attPos}.segmentedSpikeData(eList(k),goodPos{attPos}),timeRange))/diff(timeRange);
        end
        
        % ERP
        erp = squeeze(mean(lfpData{attPos}.segmentedLFPData(eList,goodPos{attPos},:),2));
        if size(erp,2)==1        
            erp=erp';       % to preserve the dimension as elecUnits x timePoints for single electrode unit case.
        end
        erpDataTMP(attPos,:,:)=erp;
        % FFT
        FFT= log10(squeeze(mean(abs(fft(lfpData{attPos}.segmentedLFPData(eList,goodPos{attPos},pos),[],3)),2)));
        if size(FFT,2)==1
            FFT=FFT';
        end
        
        fftDataTMP(attPos,:,:)=FFT;
        alphaDataTMP(attPos,:) = sum(fftDataTMP(attPos,:,alphaPos),3);
        ssvepDataTMP(attPos,:) = fftDataTMP(attPos,:,ssvepPos); %#ok<FNDSB>
    end
    
    psthData{i} = psthDataTMP;
    firingRates{i} = firingRatesTMP; 
    erpData{i} = erpDataTMP;
    fftData{i} = fftDataTMP;
    alphaData{i} = alphaDataTMP;
    ssvepData{i} = ssvepDataTMP;
    
    
end
perCorrectData=load(fullfile(folderSourceString,'Data','segmentedData','behavior',[fileNameString '_perCorrect']));  %loads the saved behavior data for each session.
perCorrect=perCorrectData.perCorrect;
end

function [colorString, colorNames] = getColorString
colorNames{1} = 'brg'; colorString{1} = 'BlueRedGreen';
colorNames{2} = 'cmy'; colorString{2} = 'CyanMagentaYellow';
end

function [fileNameStringAll,fileNameStringListArray] = getFileNameStringList

[tmpFileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

fileNameStringAll = ''; pos=1;
clear fileNameStringListArray

for i=1:length(monkeyNameList)
    for j=1:length(tmpFileNameStringList{i})
        fileNameStringAll = [cat(2,fileNameStringAll,tmpFileNameStringList{i}{j}) '|'];
        fileNameStringListArray{pos} = tmpFileNameStringList{i}(j); %#ok<*AGROW>
        pos=pos+1;
    end
end

allNames = [];
for i=1:length(monkeyNameList)
    fileNameStringAll = [cat(2,fileNameStringAll,monkeyNameList{i}) ' (N=' num2str(length(tmpFileNameStringList{i})) ')|'];
    fileNameStringListArray{pos} = tmpFileNameStringList{i};
    allNames = cat(2,allNames,tmpFileNameStringList{i});
    pos=pos+1;
end

fileNameStringAll = cat(2,fileNameStringAll,['all (N=' num2str(length(allNames)) ')']);
fileNameStringListArray{pos} = allNames;
end

function electrodeArray = getGoodElectrodes(folderSourceString,fileNameString,neuronType,populationType)

sortRatings=sortRating(fileNameString,folderSourceString);

[~,~,electrodeArrayPos]=electrodePositionOnGridMayo(1);

sua=intersect(find(sortRatings>0),find(sortRatings<4));
mua=find(sortRatings>3);
all=union(sua,mua,'sorted');

spkData{1}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H1V_StimOnset_Spikes']));
spkData{2}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H0V_StimOnset_Spikes']));
 
 if strcmp(neuronType,'SUA')
    electrodeArray{1}=intersect(sua,electrodeArrayPos(:,8:13)); % Right Array
    electrodeArray{2}=intersect(sua,electrodeArrayPos(:,1:6)); % Left Array
% electrodeArray{1} = [73 71 69 67 65 72 70 68 66 79 77 02 81 80 78 76 74 41 39 37 35 34 03 83 84 82]; % Right Array
% electrodeArray{2} = [46 50 15 17 09 11 08 90 89 55 56 57 58 52 54 19 25 12 10 92 91 94 59 60 62 21 31 29 27 20 16 22 18 96 63 61 95 32 30 28 24]; % Left Array
 
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
    plot(h,i,mData(i),'color',colorNames(i),'marker','o');
    hold(h,'on');
    errorbar(h,i,mData(i),semData(i),'color',colorNames(i));
end
set(h,'XTick',[],'XTicklabel',[]);
end

function sortRating=sortRating(fileNameString,folderSourceString)

[~,~,ratings]=xlsread(fullfile(folderSourceString,'Data','extractedData','SortRatings_CUonly'),fileNameString,'A1:A96');
if length(ratings)~=96;
    error('electrode numbers do not match')
end

for i=1:96
    
    if ischar(ratings{i})
        ratings{i}=NaN;       % to eliminate electrodes rated as x, xx and blank.
    end
end
sortRating=cell2mat(ratings);
end

function makeBehaviorPlot(h,data,colorNames)

N = size(data,3);
mData = mean(data,3);
semData = std(data,[],3)/sqrt(N);
X=(1:6);
for i=1:size(data,1)
  errorbar(h,X,mData(i,:),semData(i,:),'-o','color',colorNames(i));
  hold(h,'on'); 
end
 set(h,'XTick',1:6,'XTicklabels',1:6);
end