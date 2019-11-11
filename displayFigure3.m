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

% Figure 3a
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 3

hRaw = getPlotHandles(3,3,[0.05 0.1 0.4 0.85],0.03,0.02,0);  
hNorm = getPlotHandles(1,3,[0.5 0.1 0.45 0.85],0.03,0.02,0); linkaxes(hNorm)

% Figure 3b
set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 

hCorrWH = getPlotHandles(3,3,[0.05 0.1 0.4 0.85],0.03,0.02,0);  linkaxes(hCorrWH(1,:)); linkaxes(hCorrWH(2:3,:));
hCorrAH = getPlotHandles(3,3,[0.55 0.1 0.4 0.85],0.03,0.02,0);  linkaxes(hCorrAH)

attCueList = Figure3Dataset(1).attCueList;
attCueListCombined = Figure3Dataset(1).attCueListCombined;
nSessions = Figure3Dataset(1).nSessions;
typeList1 = [{'AOVH'} {'AOVM'} {'AIVM'} {'AIVH'} {'TOH'} {'TOM'} {'TIM'} {'TIH'} {'AOIM'} {'AOIH'} {'AIIH'} {'AIIM'}];
typeList2 = [{'H1V'} {'M1V'} {'M0V'} {'H0V'} {'H1N'} {'M1N'} {'M0N'} {'H0N'} {'M1I'} {'H1I'} {'H0I'} {'M0I'}];
cNew = [2 6 5 1 10 12 11 9 8 4 3 7];
nSessions = nSessions(cNew);
colorName = 'rbg';
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
    errorbar(hRaw(i,r),pos,mData(c,i),sData(c,i),[colorName(i),'o']);
    if i==1
        if c==1
            text(pos-0.7,mData(c,i)+sData(c,i)+0.5,['N=' num2str(nData(c,i))],'Parent',hRaw(i,r));
            text(pos-0.7,mData(c,i)-sData(c,i)-0.5,['S=' num2str(nSessions(c))],'Parent',hRaw(i,r));
        else
            text(pos-0.4,mData(c,i)+sData(c,i)+0.5,num2str(nData(c,i)),'Parent',hRaw(i,r));
            text(pos-0.4,mData(c,i)-sData(c,i)-0.5,num2str(nSessions(c)),'Parent',hRaw(i,r));
        end
    end
    hold(hRaw(i,r),'on')
    pos=pos+1;
    if mod(c,4)==0
        r = r+1; pos = 1;
    end
end
for r=1:3
    k=(r-1)*4;
    plot(hRaw(i,r),1:4,mData(k+1:k+4,i),colorName(i));
end
ylabel(hRaw(i,1),label1{i})
xlabel(hRaw(3,2),'Attention Condition')
end
YTickFRRange = [20 30];
YTickFR = 20:2:30;
YTickGammaRange = [50 150];
YTickGamma = 50:20:150;
YTickAlphaRange = [200 300];
YTickAlpha = 200:20:300;
for r=1:3
    axis(hRaw(1,r),[0 5 YTickFRRange])
    axis(hRaw(2,r),[0 5 YTickGammaRange])
    axis(hRaw(3,r),[0 5 YTickAlphaRange])
    k=(r-1)*4;
    set(hRaw(1,r),'XTick',1:4,'XTickLabel',[],'YTick',YTickFR,'YTickLabel',YTickFR,'TickDir','out','box','off');
    set(hRaw(2,r),'XTick',1:4,'XTickLabel',[],'YTick',YTickGamma,'YTickLabel',YTickGamma,'TickDir','out','box','off');
    set(hRaw(3,r),'XTick',1:4,'XTickLabel',typeList1(k+1:k+4),'XTickLabelRotation',45,'Ytick',YTickAlpha,'YTickLabel',YTickAlpha,'TickDir','out','box','off');
end

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
         errorbar(hNorm(r),pos,mData(c,i),sData(c,i),[colorName(i),'o']);
         hold(hNorm(r),'on')
         pos=pos+1;
         if mod(c,4)==0
             r = r+1; pos = 1;
         end
     end

for r=1:3
    k=(r-1)*4;
    p(i) = plot(hNorm(r),1:4,mData(k+1:k+4,i),colorName(i)); %#ok<SAGROW>
end
end
yNormRange = [-0.4 0.4];
yTickNorm = -0.4:0.1:0.4;
xNormRange = [0 5];
axis(hNorm(1),[xNormRange yNormRange])
for r=1:3
    k=(r-1)*4;
    set(hNorm(r),'XTick',1:4,'XTickLabel',typeList1(k+1:k+4),'XTickLabelRotation',45,'Ytick',yTickNorm,'YTickLabel',yTickNorm,'TickDir','out','box','off');
end
ylabel(hNorm(1),'Normalized Activity');
xlabel(hNorm(2),'Attention Condition')
legend(hNorm(1),p,'Firing Rate','Gamma Power','Alpha Power','location','northwest')
    
    
mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);
for i=1:3
x = Figure3Dataset(i).correlation_WH;
x = x(cNew);
r=1;    pos=1;
for c=1:length(x)
    mData(c,i) = mean(x{c});
    nData(c,i) = length(x{c});
    sData(c,i) = std(x{c})./sqrt(nData(c,i));
    errorbar(hCorrWH(i,r),pos,mData(c,i),sData(c,i),[colorName(i),'o']);
    if i==1
        if c==1
            text(pos-0.9,mData(c,i)+sData(c,i)+0.05,['N=' num2str(nData(c,i))],'Parent',hCorrWH(i,r));
        else
            text(pos-0.4,mData(c,i)-sData(c,i)+0.05,num2str(nData(c,i)),'Parent',hCorrWH(i,r));
        end
    end
    hold(hCorrWH(i,r),'on')
    pos=pos+1;
    if mod(c,4)==0
        r = r+1; pos = 1;
    end

end
for r=1:3
    k=(r-1)*4;
    plot(hCorrWH(i,r),1:4,mData(k+1:k+4,i),colorName(i));
end
ylabel(hCorrWH(i,1),label2{i})


end
YTickFRCRange = [0 0.4];
YTickFRC = 0:0.1:0.4;
YTickGammaCRange = [0.5 0.9];
YTickGammaC = 0.5:0.1:0.9;
YTickAlphaCRange = YTickGammaCRange;
YTickAlphaC = YTickGammaC;
for r=1:3
    k=(r-1)*4;
    axis(hCorrWH(1,1),[[0 5] YTickFRCRange]);
    axis(hCorrWH(2,1),[[0 5] YTickGammaCRange]);
    set(hCorrWH(1,r),'XTick',1:4,'XTickLabel',[],'YTick',YTickFRC,'YTickLabel',YTickFRC,'TickDir','out','box','off');
    set(hCorrWH(2,r),'XTick',1:4,'XTickLabel',[],'YTick',YTickGammaC,'YTickLabel',YTickGammaC,'TickDir','out','box','off');
    set(hCorrWH(3,r),'XTick',1:4,'XTickLabel',typeList1(k+1:k+4),'XTickLabelRotation',45,'YTick',YTickAlphaC,'YTickLabel',YTickAlphaC,'TickDir','out','box','off');
end

title(hCorrWH(1,2),'Within Hemisphere Correlation')
xlabel(hCorrWH(3,2),'Attention Condition')



mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);
for i=1:3
x = Figure3Dataset(i).correlation_AH;
x = x(cNew);
r=1; pos=1;
for c=1:length(x)
    mData(c,i) = mean(x{c});
    nData(c,i) = length(x{c});
    sData(c,i) = std(x{c})./sqrt(nData(c,i));
    errorbar(hCorrAH(i,r),pos,mData(c,i),sData(c,i),[colorName(i),'o']);
    if i==1
        if c==1
            text(pos-0.9,mData(c,i)+sData(c,i)+0.05,['N=' num2str(nData(c,i))],'Parent',hCorrAH(i,r));
        else
            text(pos-0.4,mData(c,i)-sData(c,i)+0.05,num2str(nData(c,i)),'Parent',hCorrAH(i,r));
        end
    end
    hold(hCorrAH(i,r),'on')
    pos=pos+1;
    if mod(c,4)==0
        r = r+1; pos = 1;
    end
end
for r=1:3
    k=(r-1)*4;
    plot(hCorrAH(i,r),1:4,mData(k+1:k+4,i),colorName(i));
end
ylabel(hCorrAH(i,1),label2{i})
end
xTickRange = [0 5];
yTickRange = [0 0.4];
yTick = 0:0.1:0.4;
axis(hCorrAH(1,r),[xTickRange yTickRange])
for r=1:3
    k=(r-1)*4;
    set(hCorrAH(1,r),'XTick',1:4,'XTickLabel',[],'YTick',yTick,'YTicklabel',yTick,'TickDir','out','box','off');
    set(hCorrAH(2,r),'XTick',1:4,'XTickLabel',[],'YTick',yTick,'YTicklabel',yTick,'TickDir','out','box','off');
    set(hCorrAH(3,r),'XTick',1:4,'XTickLabel',typeList2(k+1:k+4),'XTickLabelRotation',45,'YTick',yTick,'YTicklabel',yTick,'TickDir','out','box','off');
end

title(hCorrAH(1,2),'Across Hemisphere Correlation')
xlabel(hCorrAH(3,2),'Attention Condition')



