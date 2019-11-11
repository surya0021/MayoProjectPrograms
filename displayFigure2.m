function displayFigure2(SessionID,neuronType,populationType,OriChange,tapers)

close all;
if ~exist('SessionID','var');               SessionID = 'all';                      end
if ~exist('neuronType','var');              neuronType = 'All';                     end 
if ~exist('populationType','var');          populationType = 'Stimulated';          end
if ~exist('OriChange','var');               OriChange = [2 3];                      end
if ~exist('tapers','var');                  tapers = [2 3];                         end

timePeriod = 'TargetOnset';

folderSourceString = 'E:\Mayo\';
folderSave = fullfile(folderSourceString,'Data','FigureData');
makeDirectory(folderSave);
fileToSave  = fullfile(folderSave,['Figure2Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '.mat']);

if exist(fileToSave,'file')
    load(fileToSave);
else
    % Get data
    
    [ffPPC,sfPPC,freqValsPPC] ...
        = getData(SessionID,neuronType,populationType,OriChange,timePeriod,tapers);
    
    Figure2Dataset = struct('timePeriod',timePeriod,'ffPPC_WH', ffPPC{3},...
        'ffPPC_AH',ffPPC{4},'sfPPC_WH',sfPPC{3},'sfPPC_AH',sfPPC{4},...
        'freqVals',freqValsPPC);
end
    save(fileToSave,'Figure2Dataset');

% Setting Figure axes
% Figure 1
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 1

hFFPPCWH = getPlotHandles(1,2,[0.05 0.55 0.4 0.35],0.05,0.04,0);   %linkaxes(hPSTH);
hFFPPCAH = getPlotHandles(1,2,[0.55 0.55 0.4 0.35],0.05,0.04,0); %linkaxes(hERP);
hSFPPCWH= getPlotHandles(1,2,[0.05 0.1 0.4 0.35],0.05,0.04,0);   %linkaxes(hPSD);
hSFPPCAH = getPlotHandles(1,2,[0.55 0.1 0.4 0.35],0.05,0.04,0);   %linkaxes(hDelPSD);

colorNamesAttCue = 'brg';

% Plotting
for attCuePos = 1:3 % 1-Attend In, 2- Attend Out, 3- Attend Both
    
    pFFPPCWH(attCuePos) = plotData(hFFPPCWH(1),Figure2Dataset.freqVals,squeeze(Figure2Dataset.ffPPC_WH(attCuePos,:,:)),colorNamesAttCue(attCuePos));
    hold(hFFPPCWH(1),'on');
    delFFPPCWH  = squeeze(Figure2Dataset.ffPPC_WH(attCuePos,:,:)) - squeeze(Figure2Dataset.ffPPC_WH(3,:,:));
    plotData(hFFPPCWH(2),Figure2Dataset.freqVals,delFFPPCWH,colorNamesAttCue(attCuePos)); 
    hold(hFFPPCWH(2),'on');

    pFFPPCAH(attCuePos) = plotData(hFFPPCAH(1),Figure2Dataset.freqVals,squeeze(Figure2Dataset.ffPPC_AH(attCuePos,:,:)),colorNamesAttCue(attCuePos));     
    hold(hFFPPCAH(1),'on')
    delFFPPCAH  = squeeze(Figure2Dataset.ffPPC_AH(attCuePos,:,:)) - squeeze(Figure2Dataset.ffPPC_AH(3,:,:));
    plotData(hFFPPCAH(2),Figure2Dataset.freqVals,delFFPPCAH,colorNamesAttCue(attCuePos));
    hold(hFFPPCAH(2),'on') 
    
    pSFPPCWH(attCuePos) = plotData(hSFPPCWH(1),Figure2Dataset.freqVals,squeeze(Figure2Dataset.sfPPC_WH(attCuePos,:,:)),colorNamesAttCue(attCuePos));
    hold(hSFPPCWH(1),'on');
    delSFPPCWH  = squeeze(Figure2Dataset.sfPPC_WH(attCuePos,:,:)) - squeeze(Figure2Dataset.sfPPC_WH(3,:,:));
    plotData(hSFPPCWH(2),Figure2Dataset.freqVals,delSFPPCWH,colorNamesAttCue(attCuePos)); 
    hold(hSFPPCWH(2),'on') 
    
    pSFPPCAH(attCuePos) = plotData(hSFPPCAH(1),Figure2Dataset.freqVals,squeeze(Figure2Dataset.sfPPC_AH(attCuePos,:,:)),colorNamesAttCue(attCuePos)); 
    hold(hSFPPCAH(1),'on')
    delSFPPCAH  = squeeze(Figure2Dataset.sfPPC_AH(attCuePos,:,:)) - squeeze(Figure2Dataset.sfPPC_AH(3,:,:));
    plotData(hSFPPCAH(2),Figure2Dataset.freqVals,delSFPPCAH,colorNamesAttCue(attCuePos));
    hold(hSFPPCAH(2),'on') 
    
end


%%%%%%%%%%%%%%% Axis Configuration %%%%%%%%%%%%%%%%%%%%%%%555

fontSize = 10; 
tickLengthPlot = 2*get(hFFPPCWH(1),'TickLength'); 
freqRange = [0 100];
freqTicks = 0:20:100;
ffppcWHYRange = [0 1];
ffppcWHYTicks = 0:0.2:1;
delffppcWHYRange = [-0.04 0.04];
delffppcWHYTicks = -0.04:0.02:0.04;
sfppcWHYRange = [0 0.1];
sfppcWHYTicks = 0:0.02:0.1;
delsfppcWHYRange = [-0.04 0.04];
delsfppcWHYTicks = -0.04:0.02:0.04;

ffppcAHYRange = [0 0.25];
ffppcAHYTicks = 0:0.05:0.25;
delffppcAHYRange = [-0.02 0.02];
delffppcAHYTicks = -0.02:0.01:0.02;
sfppcAHYRange = [-0.005 0.025];
sfppcAHYTicks = 0:0.005:0.025;
delsfppcAHYRange = [-0.02 0.02];
delsfppcAHYTicks = -0.02:0.01:0.02;
% colorsTimePeriod = jet(3);

% FFC & delFFC Within Hemisphere Plot Graphics

% legend(hPSTH(2),'Attend inside RF','Attend outside RF','Attend Both Locations','location','best','AutoUpdate','off');

axis(hFFPPCWH(1),[freqRange ffppcWHYRange])
axis(hFFPPCWH(2),[freqRange delffppcWHYRange])

set(hFFPPCWH(1),'XTick',freqTicks,'YTick',ffppcWHYTicks,'XTicklabel',freqTicks,'YTicklabel',ffppcWHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hFFPPCWH(2),'XTick',freqTicks,'YTick',delffppcWHYTicks,'XTicklabel',freqTicks,'YTicklabel',delffppcWHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hFFPPCWH(1),{'Field-Field PPC'})
title(hFFPPCWH(2),'Change in Field-Field PPC w.r.t Neutral Condition')
ylabel(hFFPPCWH(1),'Pairwise Phase Consistency');  ylabel(hFFPPCWH(2),'Change in Pairwise Phase Consistency')
xlabel(hFFPPCWH(1),'Frequency (Hz)'); xlabel(hFFPPCWH(2),'Frequency (Hz)');

legend(hFFPPCWH(1),pFFPPCWH,'Attend inside RF','Attend outside RF','Neutral','location','southwest'); legend(hFFPPCWH(1),'boxoff');

% SFC & delSFC Within Hemisphere Plot Graphics

axis(hSFPPCWH(1),[freqRange sfppcWHYRange])
axis(hSFPPCWH(2),[freqRange delsfppcWHYRange])

set(hSFPPCWH(1),'XTick',freqTicks,'YTick',sfppcWHYTicks,'XTicklabel',freqTicks,'YTicklabel',sfppcWHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hSFPPCWH(2),'XTick',freqTicks,'YTick',delsfppcWHYTicks,'XTicklabel',freqTicks,'YTicklabel',delsfppcWHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hSFPPCWH(1),{'Spike-Field PPC'})
title(hSFPPCWH(2),'Change in Spike-Field PPC w.r.t Neutral Condition')
ylabel(hSFPPCWH(1),'Pairwise Phase Consistency');      ylabel(hSFPPCWH(2),'Change in Pairwise Phase Consistency');
xlabel(hSFPPCWH(1),'Frequency (Hz)'); xlabel(hSFPPCWH(2),'Frequency (Hz)');

legend(hSFPPCWH(1),pSFPPCWH,'Attend in RF','Attend out RF','Neutral','location','northeast'); legend(hSFPPCWH(1),'boxoff');
% FFC and delFFC Across Hemisphere Plot Graphics

axis(hFFPPCAH(1),[freqRange ffppcAHYRange])
axis(hFFPPCAH(2),[freqRange delffppcAHYRange])

set(hFFPPCAH(1),'XTick',freqTicks,'YTick',ffppcAHYTicks,'XTicklabel',freqTicks,'YTicklabel',ffppcAHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hFFPPCAH(2),'XTick',freqTicks,'YTick',delffppcAHYTicks,'XTicklabel',freqTicks,'YTicklabel',delffppcAHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hFFPPCAH(1),{'Field-Field PPC'})
title(hFFPPCAH(2),'Change in Field-Field PPC w.r.t Neutral Condition')
ylabel(hFFPPCAH(1),'Pairwise Phase Consistency');  ylabel(hFFPPCAH(2),'Change in Pairwise Phase Consistency')
xlabel(hFFPPCAH(1),'Frequency (Hz)'); xlabel(hFFPPCAH(2),'Frequency (Hz)');

legend(hFFPPCAH(2),pFFPPCAH,'Attend left','Attend right','Neutral','location','northeast'); legend(hFFPPCAH(2),'boxoff');
% SFC and delSFC Across Hemisphere Plot Graphics

axis(hSFPPCAH(1),[freqRange sfppcAHYRange])
axis(hSFPPCAH(2),[freqRange delsfppcAHYRange])

set(hSFPPCAH(1),'XTick',freqTicks,'YTick',sfppcAHYTicks,'XTicklabel',freqTicks,'YTicklabel',sfppcAHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hSFPPCAH(2),'XTick',freqTicks,'YTick',delsfppcAHYTicks,'XTicklabel',freqTicks,'YTicklabel',delsfppcAHYTicks,'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hSFPPCAH(1),{'Spike-Field PPC'})
title(hSFPPCAH(2),'Change in Spike-Field PPC w.r.t Neutral Condition')
ylabel(hSFPPCAH(1),'Pairwise Phase Consistency');      ylabel(hSFPPCAH(2),'Change in Pairwise Phase Consistency');
xlabel(hSFPPCAH(1),'Frequency (Hz)'); xlabel(hSFPPCAH(2),'Frequency (Hz)');

legend(hSFPPCAH(1),pSFPPCAH,'Attend left','Attend right','Neutral','location','northeast'); legend(hSFPPCAH(1),'boxoff');

text(60,0.8,['N = ',num2str(size(Figure2Dataset.ffPPC_WH,2))],'Color','black','FontSize',12,'parent',hFFPPCWH(1),'FontWeight','bold')
text(60,0.18,['N = ',num2str(size(Figure2Dataset.ffPPC_AH,2))],'Color','black','FontSize',12,'parent',hFFPPCAH(1),'FontWeight','bold')
text(60,0.06,['N = ',num2str(size(Figure2Dataset.sfPPC_WH,2))],'Color','black','FontSize',12,'parent',hSFPPCWH(1),'FontWeight','bold')
text(60,0.015,['N = ',num2str(size(Figure2Dataset.sfPPC_AH,2))],'Color','black','FontSize',12,'parent',hSFPPCAH(1),'FontWeight','bold')
annotation('textbox',[0.2 0.93 0.2 0.05],'String','\bf Intra Hemisphere','FontSize',12,'EdgeColor','none')
annotation('textbox',[0.7 0.93 0.2 0.05],'String','\bf Inter Hemisphere','FontSize',12,'EdgeColor','none')


annotation('textbox',[0.02 0.9 0.02 0.05],'String','\bf A','FontSize',12,'EdgeColor','none')
annotation('textbox',[0.5 0.9 0.02 0.05],'String','\bf B','FontSize',12,'EdgeColor','none')
annotation('textbox',[0.02 0.45 0.02 0.05],'String','\bf C','FontSize',12,'EdgeColor','none')
annotation('textbox',[0.5 0.45 0.02 0.05],'String','\bf D','FontSize',12,'EdgeColor','none')
end

% Separate function
function [ffPPC,sfPPC,freqValsPPC] = getData(sessionID,neuronType,populationType,OriChange,timePeriod,tapers)

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
[ffPPC,sfPPC,freqValsPPC] = loadSessionData(sessionIDList{1},neuronType,populationType,OriChange,timePeriod,tapers); % First session

if length(sessionIDList)>1
    for i = 2:numDatasets
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(sessionIDList))]);
        [ffPPCTMP,sfPPCTMP,freqValsPPCTMP] = loadSessionData(sessionIDList{i},neuronType,populationType,OriChange,timePeriod,tapers);
        
        for k=1:4 % for each array side , 3- within hemisphere combined , 4 - Across Hemisphere
            if isequal(freqValsPPC,freqValsPPCTMP)
                ffPPC{k} = cat(2,ffPPC{k},ffPPCTMP{k});
                sfPPC{k} = cat(2,sfPPC{k},sfPPCTMP{k});
            else
                error('freqVals does not match');
            end
        end
    end
end
end
function [ffPPC,sfPPC,freqValsPPC] = loadSessionData(sessionID,neuronType,populationType,OriChange,timePeriod,tapers) %#ok<STOUT>
folderSourceString = 'E:\Mayo\';
savedDatafolder = fullfile(folderSourceString,'Data','dataForFig1and2');
savedDataFile = fullfile(savedDatafolder,[sessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '_' timePeriod '_Coherence.mat']);

if exist(savedDataFile,'file')
    disp(['Loading file ' savedDataFile]);
    load(savedDataFile);
else
    error(['File not Found' savedDataFile]);
end
     
end
function p = plotData(hPlot,xs,data,colorName)

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
p = plot(hPlot,xs,mData,'color',colorName,'linewidth',1); 

end