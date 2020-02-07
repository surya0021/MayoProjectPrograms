SessionID = 'all'; 
neuronType = 'All';
populationType = 'Stimulated';
OriChange = [2 3];  
tapers = [2 3];
folderSourceString = 'E:/Mayo';
folderSave = fullfile(folderSourceString,'Data','FigureData');
fileToSave  = fullfile(folderSave,['Figure3Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '.mat']);
dataType{1} = 'Firing Rate';
dataType{2} = 'Gamma Power';
dataType{3} = 'Alpha Power';

if exist(fileToSave,'file')
    load(fileToSave)
else
    for dataTypeNum = 1:3
        disp(['Working on ' dataType{dataTypeNum}]);
        mCombinedData = cell(1,2);  rData = cell(1,2);  rCombinedData = cell(1,2);
        for processType = 1:2
            [~,mCombinedData{processType},rData{processType},rCombinedData{processType},typeList1,typeList2,nSessions] = displayHvsMAnalysisMagnitudesCorrelations2(dataTypeNum,processType);
        end
        Figure3Dataset(dataTypeNum) = struct('dataType',dataType{dataTypeNum},'rawData', mCombinedData(1),...
            'zScoreData',mCombinedData(2),'correlation_WH',rCombinedData(1),'correlation_AH',{rData{1}(:,3)'},...
            'attCueList',{typeList1},'attCueListCombined',{typeList2},'nSessions',nSessions); %#ok<SAGROW>
    end
    
    save(fileToSave,'Figure3Dataset')
end

% Figure 3
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 3

hRaw = getPlotHandles(3,1,[0.05 0.1 0.4 0.85],0.03,0.02,0);  
hNorm = getPlotHandles(1,1,[0.5 0.1 0.45 0.85],0.03,0.02,0); 


attCueList = Figure3Dataset(1).attCueList;
attCueListCombined = Figure3Dataset(1).attCueListCombined;
nSessions = Figure3Dataset(1).nSessions;
typeList1 = [{'AOVH'} {'AOVM'} {'AIVM'} {'AIVH'} {'TOH'} {'TOM'} {'TIM'} {'TIH'} {'AOIM'} {'AOIH'} {'AIIH'} {'AIIM'}];
typeList2 = [{'H1V'} {'M1V'} {'M0V'} {'H0V'} {'H1N'} {'M1N'} {'M0N'} {'H0N'} {'M1I'} {'H1I'} {'H0I'} {'M0I'}];
cNew = [2 6 5 1 10 12 11 9 8 4 3 7];
xIndex = [1:4 6:9 11:14];
nSessions = nSessions(cNew);
colorName = gray(4);
label1 = [{'Raw Firing Rate (spikes/s)'} {'Raw Gamma Power (\muV^2)'} {' Raw Alpha Power (\muV^2)'}];
label2 = [{'Firing Rate'} {'Gamma Power'} {'Alpha Power'}];
%Plot
% Raw Data
mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);
for i=1:3
x = Figure3Dataset(i).rawData;
x = x(cNew);
r = 1;  pos = 1;
for c=1:length(x)
    mData(c,i) = mean(x{c});
    nData(c,i) = length(x{c});
    sData(c,i) = std(x{c})./sqrt(nData(c,i));
    errorbar(hRaw(i),pos,mData(c,i),sData(c,i),'Color',colorName(i,:),'Marker','o');
    if i==1
        if c==1
            text(pos-0.7,mData(c,i)+sData(c,i)+0.5,['N=' num2str(nData(c,i))],'Parent',hRaw(i));
            text(pos-0.7,mData(c,i)-sData(c,i)-0.5,['S=' num2str(nSessions(c))],'Parent',hRaw(i));
        else
            text(pos-0.4,mData(c,i)+sData(c,i)+0.5,num2str(nData(c,i)),'Parent',hRaw(i));
            text(pos-0.4,mData(c,i)-sData(c,i)-0.5,num2str(nSessions(c)),'Parent',hRaw(i));
        end
    end
    hold(hRaw(i),'on')
    if mod(c,4)==0
        pos = pos+2;
    else
        pos = pos+1;
    end
end
for j=1:3
    r = (j-1)*4;
    plot(hRaw(i),xIndex(r+1:r+4),mData(r+1:r+4,i),'Color',colorName(i,:),'LineWidth',1.2);
end

line([0 15],[mData(4,i) mData(4,i)],'Color',colorName(i,:),'Parent',hRaw(i)) 
line([0 15],[mData(1,i) mData(1,i)],'Color',colorName(i,:),'Parent',hRaw(i)) 

ylabel(hRaw(i),label1{i})
xlabel(hRaw(3),'Attention Condition')
end
YTickFRRange = [18 26];
YTickFR = 18:2:26;
YTickGammaRange = [70 110];
YTickGamma = 70:10:110;
YTickAlphaRange = [200 300];
YTickAlpha = 200:20:300;

axis(hRaw(1),[0 15 YTickFRRange])
axis(hRaw(2),[0 15 YTickGammaRange])
axis(hRaw(3),[0 15 YTickAlphaRange])

set(hRaw(1),'XTick',xIndex,'XTickLabel',[],'YTick',YTickFR,'YTickLabel',YTickFR,'TickDir','out','box','off');
set(hRaw(2),'XTick',xIndex,'XTickLabel',[],'YTick',YTickGamma,'YTickLabel',YTickGamma,'TickDir','out','box','off');
set(hRaw(3),'XTick',xIndex,'XTickLabel',typeList1,'XTickLabelRotation',45,'Ytick',YTickAlpha,'YTickLabel',YTickAlpha,'TickDir','out','box','off');


% NormData
 mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);
 for i=1:3
     x = Figure3Dataset(i).zScoreData;
     x = x(cNew);
     r=1; pos=1;
     for c=1:length(x)
         mData(c,i) = mean(x{c});
         nData(c,i) = length(x{c});
         sData(c,i) = std(x{c})./sqrt(nData(c,i));
         errorbar(hNorm,pos,mData(c,i),sData(c,i),'Color',colorName(i,:),'Marker','o');
         hold(hNorm,'on')
         if mod(c,4)==0
             pos = pos+2;
         else
             pos = pos+1;
         end
     end
     
     for j=1:3
         r=(j-1)*4;
         p(i) = plot(hNorm,xIndex(r+1:r+4),mData(r+1:r+4,i),'Color',colorName(i,:),'LineWidth',1.2); %#ok<SAGROW>
     end
     line([0 15],[mData(4,i) mData(4,i)],'Color',colorName(i,:),'Parent',hNorm)
 end

line([0 15],[mData(1,1) mData(1,i)],'Color',colorName(1,:),'LineStyle','--','Parent',hNorm)
yNormRange = [-0.2 0.4];
yTickNorm = -0.2:0.1:0.4;
xNormRange = [0 15];
axis(hNorm,[xNormRange yNormRange])
set(hNorm,'XTick',xIndex,'XTickLabel',typeList1,'XTickLabelRotation',45,'Ytick',yTickNorm,'YTickLabel',yTickNorm,'TickDir','out','box','off');
ylabel(hNorm,'Normalized Activity');
xlabel(hNorm,'Attention Condition')
legend(hNorm,p,'Firing Rate','Gamma Power','Alpha Power','location','northwest')
    
annotation('textbox',[0.02 0.95 0.02 0.05],'String','\bf A','FontSize',12,'EdgeColor','none')
annotation('textbox',[0.47 0.95 0.02 0.05],'String','\bf B','FontSize',12,'EdgeColor','none')    



