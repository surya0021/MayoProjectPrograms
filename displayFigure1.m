function displayFigure1(SessionID,neuronType,populationType,OriChange,tapers)

close all;
if ~exist('SessionID','var');               SessionID = 'all';                      end
if ~exist('neuronType','var');              neuronType = 'All';                     end 
if ~exist('populationType','var');          populationType = 'Stimulated';          end
if ~exist('OriChange','var');               OriChange = [2 3];                      end
if ~exist('tapers','var');                  tapers = [2 3];                         end
timePeriod{1} = 'Baseline';
timePeriod{2} = 'StimOnset';
timePeriod{3} = 'TargetOnset';

folderSourceString = 'E:\Mayo\';
folderSave = fullfile(folderSourceString,'Data','FigureData');
makeDirectory(folderSave);
fileToSave  = fullfile(folderSave,['Figure1Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '.mat']);

if exist(fileToSave,'file')
    load(fileToSave);
else
    % Get data combined across Hemispheres
    for iTimePeriod = 1:length(timePeriod)
        clear psthData xsFR erpData timeVals psdData freqVals
        
        [psthData,xsFR,erpData,timeVals,psdData,freqVals] ...
           = getData(SessionID,neuronType,populationType,OriChange,timePeriod{iTimePeriod},tapers);

        Figure1Dataset(iTimePeriod) = struct('timePeriod',timePeriod{iTimePeriod},'psthData', psthData{3},...
                                            'xsFR',xsFR,'erpData',erpData{3},'timeVals',timeVals,...
                                            'psdData',psdData{3},'freqVals',freqVals);
    end
    save(fileToSave,'Figure1Dataset');
end

% Setting Figure axes
% Figure 1
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 1

hPSTH = getPlotHandles(1,2,[0.05 0.72 0.4 0.25],0.01,0.02,0);   %linkaxes(hPSTH);
hERP = getPlotHandles(1,2,[0.55 0.72 0.4 0.25],0.01,0.02,0); %linkaxes(hERP);
hPSD = getPlotHandles(1,3,[0.05 0.06 0.4 0.55],0.05,0.02,0);   linkaxes(hPSD);
hDelPSD = getPlotHandles(3,1,[0.55 0.06 0.4 0.55],0.02,0.02,0);   linkaxes(hDelPSD);

colorNamesAttCue = 'brg';

% Plotting
for attCuePos = 1:3 % 1-Attend In, 2- Attend Out, 3- Attend Both
    pPSTH(attCuePos) = plotData(hPSTH(1),Figure1Dataset(2).xsFR,squeeze(Figure1Dataset(2).psthData(attCuePos,:,:)),colorNamesAttCue(attCuePos)); % StimOnset
    hold(hPSTH(1),'on');
    plotData(hPSTH(2),Figure1Dataset(3).xsFR,squeeze(Figure1Dataset(3).psthData(attCuePos,:,:)),colorNamesAttCue(attCuePos)); % TargetOnset
    hold(hPSTH(2),'on');

    pERP(attCuePos) = plotData(hERP(1),Figure1Dataset(2).timeVals,squeeze(Figure1Dataset(2).erpData(attCuePos,:,:)),colorNamesAttCue(attCuePos)); % StimOnset
    hold(hERP(1),'on')    
    plotData(hERP(2),Figure1Dataset(3).timeVals,squeeze(Figure1Dataset(3).erpData(attCuePos,:,:)),colorNamesAttCue(attCuePos)); % TargetOnset
    hold(hERP(2),'on')    
end

for i = 1:3
    for attCuePos = 1:3
        PSD = log10(squeeze(Figure1Dataset(i).psdData(attCuePos,:,:)));
        delPSD  = 10*(log10(squeeze(Figure1Dataset(i).psdData(attCuePos,:,:))) - log10(squeeze(Figure1Dataset(i).psdData(3,:,:))));
        
        pPSD{i}(attCuePos) = plotData(hPSD(i),Figure1Dataset(i).freqVals,PSD,colorNamesAttCue(attCuePos));
        hold(hPSD(i),'on');
        
        pDelPSD{i}(attCuePos) = plotData(hDelPSD(i),Figure1Dataset(i).freqVals,delPSD,colorNamesAttCue(attCuePos));
        hold(hDelPSD(i),'on');
    end
end

%%%%%%%%%%%%%%% Axis Configuration %%%%%%%%%%%%%%%%%%%%%%%555

fontSize = 10;
FRTicks = 0:20:40; 
timeTicks{1} = [-0.2 0 0.2 0.4];
timeTicks{2} = [-0.4 -0.2 0];
tickLengthPlot = 2*get(hPSTH(1),'TickLength'); 
timeVals{1} = [-0.25 0.5]; timeVals{2} = [-0.5 0.1];
frRange = [0 50];
erpRange = [-300 50];
erpTicks = [-250 -150 0];

% colorsTimePeriod = jet(3);

% PSTH & ERP Plot Graphics
text(0.25,40,['N = ',num2str(size(Figure1Dataset(1).psthData,2))],'Color','black','FontSize',14,'parent',hPSTH(1),'FontWeight','bold')
legend(hPSTH(2),pPSTH,'Attend inside RF','Attend outside RF','Neutral','location','best'); legend(hPSTH(2),'boxoff');
legend(hERP(2),pERP,'Attend inside RF','Attend outside RF','Neutral','location','southwest'); legend(hERP(2),'boxoff');
for i = 1:2
    axis(hPSTH(i),[timeVals{i} frRange]); 
    axis(hERP(i),[timeVals{i} erpRange])
    line([0 0],[0 50],'color','black','lineWidth',2,'parent',hPSTH(i));
    line([0 0],[-300 50],'color','black','lineWidth',2,'parent',hERP(i));
end
set(hPSTH(1),'XTick',timeTicks{1},'YTick',FRTicks,'YTicklabel',FRTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPSTH(2),'XTick',timeTicks{2},'YTick',FRTicks,'YTicklabel',FRTicks,'YAxisLocation','right','fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
ylabel(hPSTH(1),{'\bf Stimulus Onset'; '\rm Firing Rate (spikes/s)'})
ylabel(hPSTH(2),'\bf Target Onset')
xlabel(hPSTH(1),'Time w.r.t StimOnset (s)'); xlabel(hPSTH(2),'Time w.r.t TargetOnset (s)');

set(hERP(1),'XTick',timeTicks{1},'YTick',erpTicks,'YTicklabel',erpTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hERP(2),'XTick',timeTicks{2},'YTick',erpTicks,'YTicklabel',erpTicks,'YAxisLocation','right','fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
ylabel(hERP(1),{'\bf Stimulus Onset'; '\rm ERP (\muV)'});
ylabel(hERP(2),'\bf TargetOnset');
xlabel(hERP(1),'Time w.r.t Stimulus Onset (s)');xlabel(hERP(2),'Time w.r.t Target Onset (s)');


patchLocX{1} = [-0.25 0 0 -0.25]; 
patchLocX{2} = [0.25 0.5 0.5 0.25];
patchLocX{3} = [-0.5 0 0 -0.5]; 
patchLocY{1} = [0 0 1 1];
patchLocY{2} = [-300 -300 -295 -295];
for j=1:2
patch(patchLocX{j},patchLocY{1},'k','parent',hPSTH(1))
patch(patchLocX{j},patchLocY{2},'k','parent',hERP(1))
end
patch(patchLocX{3},patchLocY{1},'k','parent',hPSTH(2))
patch(patchLocX{3},patchLocY{2},'k','parent',hERP(2))

% PSD and delPSD Plot Graphics

freqRange = [0 100];
psdYRange = [-0.2 4];
delPSDYRange = [-0.6 0.6];
freqTicks = 0:20:100;
psdYTicks = 0:1:4;
delPSDYTicks = -0.6:0.3:0.6;

axis(hPSD(1),[freqRange psdYRange]);    
axis(hDelPSD(1),[freqRange delPSDYRange]);
for i=1:2
set(hPSD(i),'XTick',freqTicks,'YTick',psdYTicks,'XTicklabel',[],'YTicklabel',psdYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hDelPSD(i),'XTick',freqTicks,'YTick',delPSDYTicks,'XTicklabel',[],'YTicklabel',delPSDYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
set(hPSD(3),'XTick',freqTicks,'YTick',psdYTicks,'XTicklabel',freqTicks,'YTicklabel',psdYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hDelPSD(3),'XTick',freqTicks,'YTick',delPSDYTicks,'XTicklabel',freqTicks,'YTicklabel',delPSDYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')

set(hPSTH(2),'XTick',timeTicks{2},'YTick',FRTicks,'YTicklabel',FRTicks,'YAxisLocation','right','fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')

xlabel(hPSD(3),'Frequency (Hz)');
xlabel(hDelPSD(3),'Frequency (Hz)');
tpStr{1} = 'Baseline';   tpStr{2} = 'Stimulus Onset';   tpStr{3} = 'Target Onset';
for i=1:3
    ylabel(hPSD(i),{['\bf' tpStr{i}] ;'\rm Power log_{10}(\muV^2)'});
    ylabel(hDelPSD(i),{['\bf' tpStr{i}] ;'\rm Change in Power (dB)'});
end
%title
title(hPSD(2),'Power Spectral Density(PSD)');
title(hDelPSD(1), 'Change in PSD w.r.t Neutral Condition')

%legend
legend(hPSD(1),pPSD{1},'Attend inside RF','Attend outside RF','Neutral','location','northeast');    legend(hPSD(1),'boxoff');
legend(hDelPSD(1),pDelPSD{1},'Attend inside RF','Attend outside RF','Neutral','location','northeast');  legend(hDelPSD(1),'boxoff');

end

% Separate function
function [psthData,xsFR,erpData,timeVals,psdData,freqVals] = getData(sessionID,neuronType,populationType,OriChange,timePeriod,tapers)

sessionNameString = getAttentionExperimentDetails;
if strcmp(sessionID,'all')
    sessionIDList = cat(2,sessionNameString{1},sessionNameString{2});
elseif strcmp(sessionID,'arturo')
    sessionIDList = sessionNameString{1};
elseif stcmp(sessionID,'wiggin')
    sessionIDList = sessionNameString{2};
end

numDatasets = length(sessionIDList);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[psthData,xsFR,erpData,timeVals,psdData,freqVals] = loadSessionData(sessionIDList{1},neuronType,populationType,OriChange,timePeriod,tapers); % First session

if length(sessionIDList)>1
    for i = 2:numDatasets
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(sessionIDList))]);
        [psthDataTMP,xsFRTMP,erpDataTMP,timeValsTMP,psdDataTMP,freqValsTMP] = loadSessionData(sessionIDList{i},neuronType,populationType,OriChange,timePeriod,tapers);
        
        for k=1:3 % for each array side and 3- both arrays combined
            if isequal(xsFR,xsFRTMP)
                psthData{k} = cat(2,psthData{k},psthDataTMP{k});
            else
                error('xsFR does not match');
            end
            if isequal(timeVals,timeValsTMP)
                erpData{k} = cat(2,erpData{k},erpDataTMP{k});
            else
                error('timeVals do not match');
            end
            if isequal(freqVals,freqValsTMP)
                psdData{k} = cat(2,psdData{k},psdDataTMP{k});
            else
                error('freqVals do not match');
            end
        end
    end
end
end
function [psthData,xsFR,erpData,timeVals,psdData,freqVals] = loadSessionData(sessionID,neuronType,populationType,OriChange,timePeriod,tapers) %#ok<STOUT>
folderSourceString = 'E:\Mayo\';
savedDatafolder = fullfile(folderSourceString,'Data','dataForFig1and2');
savedDataFile = fullfile(savedDatafolder,[sessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '_' timePeriod '.mat']);

if exist(savedDataFile,'file')
    disp(['Loading file ' savedDataFile]);
    load(savedDataFile);
else
    error(['File not Found' savedDataFile]);
end
     
end
function p= plotData(hPlot,xs,data,colorName)

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

patch(xsLong,ysLong,colorName2,'FaceAlpha',0.2,'parent',hPlot);
hold(hPlot,'on');
p=plot(hPlot,xs,mData,'color',colorName,'linewidth',1); 

end