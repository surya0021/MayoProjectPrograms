function displayFigure1_Trial(SessionID,neuronType,populationType,trialOutcome,OrientationChange)

close all;
if ~exist('SessionID','var');               SessionID = 'all';                      end
if ~exist('neuronType','var');              neuronType = 'All';                     end 
if ~exist('populationType','var');          populationType = 'Stimulated';          end
if ~exist('trialOutcome','var');            trialOutcome = 'H';                     end
if ~exist('OrientationChange','var');       OrientationChange = 'Selected';         end

timePeriod{1} = 'Baseline';
timePeriod{2} = 'StimOnset';
timePeriod{3} = 'TargetOnset';

folderSourceString = 'E:\Projects\MayoProject\';
folderSave = fullfile(folderSourceString,'Data','FigureData');
makeDirectory(folderSave);
fileToSave  = fullfile(folderSave,['Figure1Dataset_' SessionID neuronType populationType trialOutcome OrientationChange '.mat']);


if exist(fileToSave,'file')
    load(fileToSave);
else
    % Get data combined across Hemispheres
    for iTimePeriod = 1:length(timePeriod)
        clear psthData xsFR erpData timeVals fftData
        clear firingRates alphaData SSVEPData

        [psthData,xsFR,erpData,timeVals,fftData,freqVals,firingRates,alphaData,ssvepData] ...
           = getData(SessionID,neuronType,populationType,trialOutcome,OrientationChange,timePeriod{iTimePeriod});

        Figure1Dataset(iTimePeriod) = struct('timePeriod',timePeriod{iTimePeriod},'psthData', psthData{3}([1,5,2],:,:),...
                                            'xsFR',xsFR,'erpData',erpData{3}([1,5,2],:,:),'timeVals',timeVals,...
                                            'fftData',fftData{3}([1,5,2],:,:),'freqVals',freqVals,...
                                            'firingRates',firingRates{3}([1,5,2],:,:),'alphaData',alphaData{3}([1,5,2],:,:),...
                                            'ssvepData',ssvepData{3}([1,5,2],:,:)); % Data Rearranged to keep 3 Attention conditions- AttendIn-AttendBoth-AttendOut
    end
    save(fileToSave,'Figure1Dataset');
end

% Setting Figure axes
% Figure 1
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 1

hPSTH = getPlotHandles(1,2,[0.075 0.72 0.4 0.25],0.01,0.02,0);%linkaxes(hPSTH);
hERP = getPlotHandles(1,2,[0.075 0.45 0.4 0.25],0.01,0.02,0);%linkaxes(hERP);
hFFT = getPlotHandles(1,3,[0.075 0.08 0.4 0.25],0.01,0.02,0);linkaxes(hFFT);

hBarFR = getPlotHandles(1,3,[0.55 0.72 0.4 0.25],0.01,0.02,0);linkaxes(hBarFR);
hBarAlpha = getPlotHandles(1,3,[0.55 0.45 0.4 0.25],0.01,0.02,0);linkaxes(hBarAlpha)
hBarSSVEP = getPlotHandles(1,3,[0.55 0.18 0.4 0.25],0.01,0.02,0);linkaxes(hBarSSVEP);

colorNamesAttCue = 'bgr';

% Plotting
for attCuePos = 1:3 % 1-Attend In, 2- Attend Both, 3- Attend Out
    plot(hPSTH(1),Figure1Dataset(2).xsFR,squeeze(mean(Figure1Dataset(2).psthData(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos)); % StimOnset
    hold(hPSTH(1),'on');
    plot(hPSTH(2),Figure1Dataset(3).xsFR,squeeze(mean(Figure1Dataset(3).psthData(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos)); % TargetOnset
    hold(hPSTH(2),'on');

    plot(hERP(1),Figure1Dataset(2).timeVals,squeeze(mean(Figure1Dataset(2).erpData(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos)); % StimOnset
    hold(hERP(1),'on')    
    plot(hERP(2),Figure1Dataset(3).timeVals,squeeze(mean(Figure1Dataset(3).erpData(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos)); % TargetOnset
    hold(hERP(2),'on')    
end

for i = 1:3
    for attCuePos = 1:3
        plot(hFFT(i),Figure1Dataset(i).freqVals,squeeze(mean(Figure1Dataset(i).fftData(attCuePos,:,:),2)),'color',colorNamesAttCue(attCuePos));
        hold(hFFT(i),'on');
    end
    makeBarPlot(hBarFR(i),Figure1Dataset(i).firingRates,colorNamesAttCue);
    makeBarPlot(hBarAlpha(i),Figure1Dataset(i).alphaData,colorNamesAttCue);
    makeBarPlot(hBarSSVEP(i),Figure1Dataset(i).ssvepData,colorNamesAttCue);
end

%%%%%%%%%%%%%%% Axis Configuration %%%%%%%%%%%%%%%%%%%%%%%555

fontSize = 12;
FRTicks = 0:20:40; 
timeTicks{1} = [-0.2 0 0.2 0.4];
timeTicks{2} = [-0.4 -0.2 0];
tickLengthPlot = 2*get(hPSTH(1),'TickLength'); 
timeVals{1} = [-0.25 0.5]; timeVals{2} = [-0.5 0.1];
frRange = [0 50];
erpRange = [-300 50];
erpTicks = [-250 -150 0];
freqTicks = 10:10:50;
fftYTicks = 3.5:0.5:4.5;
freqRange = [0 50];
fftYRange = [3 5];
% colorsTimePeriod = jet(3);

% PSTH & ERP Plot Graphics
text(0.25,40,['N = ',num2str(size(Figure1Dataset(1).alphaData,2))],'Color','black','FontSize',14,'parent',hPSTH(1),'FontWeight','bold')
legend(hPSTH(2),'Attend inside RF','Attend outside RF','Attend Both Locations','location','best','AutoUpdate','off');
for i = 1:2
    axis(hPSTH(i),[timeVals{i} frRange])
    axis(hERP(i),[timeVals{i} erpRange])
    line([0 0],[0 50],'color','black','lineWidth',2,'parent',hPSTH(i));
    line([0 0],[-300 50],'color','black','lineWidth',2,'parent',hERP(i));
end
set(hPSTH(1),'XTick',timeTicks{1},'YTick',FRTicks,'XTicklabel',[],'YTicklabel',FRTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPSTH(2),'XTick',timeTicks{2},'YTick',FRTicks,'XTicklabel',[],'YTicklabel',[],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hPSTH(1),'Stimulus Onset'); title(hPSTH(2),'Target Onset')
ylabel(hPSTH(1),'Firing Rate (spikes/s)')

set(hERP(1),'XTick',timeTicks{1},'YTick',erpTicks,'YTicklabel',erpTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hERP(2),'XTick',timeTicks{2},'YTick',erpTicks,'YTicklabel',[],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
ylabel(hERP(1),'ERP (\muV)');
xlabel(hERP(1),'Time w.r.t StimOnset (s)');xlabel(hERP(2),'Time w.r.t TargetOnset (s)');


patchLocX{1} = [-0.25 0 0 -0.25]; 
patchLocX{2} = [0.25 0.5 0.5 0.25];
patchLocY{1} = [0 0 2 2];
patchLocY{2} = [-300 -300 -286 -286];
for j=1:2
patch(patchLocX{j},patchLocY{1},'k','parent',hPSTH(1))
patch(patchLocX{j},patchLocY{2},'k','parent',hERP(1))
end
patch(patchLocX{1},patchLocY{1},'k','parent',hPSTH(2))
patch(patchLocX{1},patchLocY{2},'k','parent',hERP(2))

% FFT and Bar Plot Graphics
axis(hFFT(1),[freqRange fftYRange]);
ylabel(hFFT(1),'Log FFT');
set(hBarFR(1),'YLim',[0 30]);
set(hBarAlpha(1),'YLim',[7.5 8.5]);
set(hBarSSVEP(1),'YLim',[3.7 4]);
ylabel(hBarFR(1),'Firing Rate (spikes/s)'); 
ylabel(hBarAlpha(1),'Alpha, log_1_0(FFT)');
ylabel(hBarSSVEP(1),'SSVEP, log_1_0(SSVEP)');

for i = 1:3
    title(hFFT(i),timePeriod{i});
    xlabel(hFFT(i),'Frequency (Hz)')
    title(hBarFR(i),timePeriod{i});
    if i==1
        set(hFFT(i),'XTick',freqTicks,'YTick',fftYTicks,'YTicklabel',fftYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hBarFR(i),'XTick',[],'YTick',[5 15 25],'YTicklabel',[5 15 25],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hBarAlpha(i),'XTick',[],'YTick',[7.8 8.1 8.4],'YTicklabel',[7.8 8.1 8.4],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hBarSSVEP(i),'XTick',[],'YTick',[3.75 3.85 3.95],'YTicklabel',[3.75 3.85 3.95],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    else
        set(hFFT(i),'XTick',freqTicks,'YTick',fftYTicks,'YTicklabel',[],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hBarFR(i),'XTick',[],'YTick',[5 15 25],'YTicklabel',[],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hBarAlpha(i),'XTick',[],'YTick',[7.4 7.9 8.4],'YTicklabel',[],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hBarSSVEP(i),'XTick',[],'YTick',[3.75 3.85 3.95],'YTicklabel',[],'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    end
end
end

% Separate function
function [psthData,xsFR,erpData,timeVals,fftData,freqVals,firingRates,alphaData,ssvepData] = getData(SessionID,neuronType,populationType,trialOutcome,OrientationChange,timePeriod)
                
        [fileNameStringList,monkeyNameList] = getAttentionExperimentDetails;
        if strcmp(SessionID,'all')
            pos = 1;
            for i=1:length(monkeyNameList)
                for j=1:length(fileNameStringList{i})
                    SessionIDList{pos} = fileNameStringList{i}(j); %#ok<*AGROW>
                    pos=pos+1;
                end
            end
        elseif strcmp(SessionID,'arturo')
            SessionIDList = fileNameStringList{1};
        elseif strcmp(SessionID,'wiggin')
            SessionIDList = fileNameStringList{2};
        end
        
        numDatasets = length(SessionIDList);
        disp(['Working on dataset 1 of ' num2str(numDatasets)]);
        [psthData,xsFR,erpData,timeVals,fftData,freqVals,firingRates,alphaData,ssvepData] = loadSessionData(SessionIDList{1},neuronType,populationType,trialOutcome,OrientationChange,timePeriod); % First session

if length(SessionIDList)>1
for i = 2:numDatasets
    disp(['Working on dataset ' num2str(i) ' of ' num2str(length(SessionIDList))]);
    [psthDataTMP,xsFRTMP,erpDataTMP,timeValsTMP,fftDataTMP,freqValsTMP,firingRatesTMP,alphaDataTMP,ssvepDataTMP] = loadSessionData(SessionIDList{i},neuronType,populationType,trialOutcome,OrientationChange,timePeriod); 
    
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
end
end
end
function [psthData,xsFR,erpData,timeVals,fftData,freqVals,firingRates,alphaData,ssvepData] = loadSessionData(SessionID,neuronType,populationType,trialOutcome,OrientationChange,timePeriod) %#ok<STOUT>
folderSourceString = 'E:\Projects\MayoProject\';
savedDatafolder = fullfile(folderSourceString,'Data','savedData');
savedDataFile = fullfile(savedDatafolder,[SessionID{1} neuronType populationType trialOutcome OrientationChange timePeriod '.mat']);

if exist(savedDataFile,'file')
    disp(['Loading file ' savedDataFile]);
    load(savedDataFile);
else
    error(['File not Found' savedDataFile]);
end
     
end
function makeBarPlot(h,data,colorNames)

N = size(data,2);
mData = mean(data,2);
semData = std(data,[],2)/sqrt(N);

for i=1:size(data,1)
    bar(h,i,mData(i),colorNames(i),'lineWidth',1.5);
    hold(h,'on');
    errorbar(h,i,mData(i),semData(i),'color','k','marker','.','lineStyle','-','lineWidth',1.5);
end
set(h,'XTick',[],'XTicklabel',[]);
end