SessionID = 'all'; 
neuronType = 'All';
populationType = 'Stimulated';
OriChange = [2 3];  
tapers = [2 3];
combineFlag = 1;
% trialCutoff = 40;
folderSourceString = 'E:/Mayo';
folderSave = fullfile(folderSourceString,'Data','FigureData');
fileToSave  = fullfile(folderSave,['Figure3Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '.mat']);
% fileToSave  = fullfile(folderSave,['Figure3Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '_trialCutoff' num2str(trialCutoff) '.mat']);
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
            [~,mCombinedData{processType},rData{processType},rCombinedData{processType},typeList1,typeList2,nSessions] = displayHvsMAnalysisMagnitudesCorrelations2(dataTypeNum,processType,trialCutoff);
        end
        Figure3Dataset(dataTypeNum) = struct('dataType',dataType{dataTypeNum},'rawData', mCombinedData(1),...
            'zScoreData',mCombinedData(2),'correlation_WH',rCombinedData(1),'correlation_AH',{rData{1}(:,3)'},...
            'attCueList',{typeList1},'attCueListCombined',{typeList2},'nSessions',nSessions); %#ok<SAGROW>
    end
    
    save(fileToSave,'Figure3Dataset')
end

set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 

hCorrWH = getPlotHandles(1,1,[0.05 0.1 0.4 0.85],0.03,0.02,0);  linkaxes(hCorrWH)
hCorrAH = getPlotHandles(1,1,[0.55 0.1 0.4 0.85],0.03,0.02,0);  linkaxes(hCorrAH)

attCueList = Figure3Dataset(1).attCueList;
nSessions = Figure3Dataset(1).nSessions;
attCueListCombined = Figure3Dataset(1).attCueListCombined;
typeList1 = [{'AOVH'} {'AOVM'} {'AIVM'} {'AIVH'} {'TOH'} {'TOM'} {'TIM'} {'TIH'} {'AOIM'} {'AOIH'} {'AIIH'} {'AIIM'}];
xIndex1 = [1:4 6:9 11:14];

if combineFlag
    typeList2 = [{'HV'} {'MV'} {'HN'} {'MN'} {'MI'} {'HI'}];
    xIndex2 = [1 2 4 5 7 8];
else
    typeList2 = [{'H1V'} {'M1V'} {'M0V'} {'H0V'} {'H1N'} {'M1N'} {'M0N'} {'H0N'} {'M1I'} {'H1I'} {'H0I'} {'M0I'}];
    xIndex2 = xIndex1;
end
cNew = [2 6 5 1 10 12 11 9 8 4 3 7];

nSessions = nSessions(cNew);
colorName = gray(4);

% Intra-Hemispheric Correlation
mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);  

for i=1:3
    x = Figure3Dataset(i).correlation_WH;
    x = x(cNew);
    pos=1;
    for c=1:12
        z{c} = x{c}-mean(x{1}); %#ok<SAGROW>
        mData(c,i) = mean(z{c});
        nData(c,i) = length(z{c});
        sData(c,i) = std(z{c})./sqrt(nData(c,i));
        errorbar(hCorrWH,pos,mData(c,i),sData(c,i),'Color',colorName(i,:),'Marker','o');
        if i==2
            if c==1
                text(pos-0.9,mData(c,i)+sData(c,i)+0.005,['N=' num2str(nData(c,i))],'Parent',hCorrWH);
                text(pos-0.9,mData(c,i)-sData(c,i)-0.005,['S=' num2str(nSessions(c))],'Parent',hCorrWH);
            else
                text(pos-0.4,mData(c,i)+sData(c,i)+0.005,num2str(nData(c,i)),'Parent',hCorrWH);
                text(pos-0.2,mData(c,i)-sData(c,i)-0.005,num2str(nSessions(c)),'Parent',hCorrWH);
            end
        end
        hold(hCorrWH,'on')
        
        if mod(c,4)==0
            pos = pos+2;
        else
            pos = pos+1;
        end
    end
    for r=1:3
        k=(r-1)*4;
        p(i) = plot(hCorrWH,xIndex1(k+1:k+4),mData(k+1:k+4,i),'Color',colorName(i,:),'LineWidth',1.2);
    end
    
end
YTickCWHRange = [-0.08 0.08];
YTickCWH = -0.08:0.02:0.08;
axis(hCorrWH,[[0 15] YTickCWHRange]);
set(hCorrWH,'XTick',xIndex1,'XTickLabel',typeList1,'YTick',YTickCWH,'YTickLabel',YTickCWH,'XTickLabelRotation',45,'TickDir','out','box','off');

title(hCorrWH,'Intra-Hemisphere')
xlabel(hCorrWH,'Attention Condition')
ylabel(hCorrWH,'Change in correlation')
legend(hCorrWH,p,'Firing Rate','Gamma Power','Alpha Power','location','northeast')

% Inter-Hemispheric Correlation
if combineFlag
    numCondition = 6;
    grp = 2;
else
    numCondition = 12;
    grp = 4;
end
mData = zeros(numCondition,3); sData = zeros(numCondition,3);   nData = zeros(numCondition,3);
clear x
for i=1:3
    corrAH = Figure3Dataset(i).correlation_AH;
    if combineFlag
        m = 1;
        for j=1:6
            x{j}= cat(1,corrAH{m},corrAH{m+1});
            m = m+2;
        end
        order = [1 3 5 6 4 2];
        x = x(order);
    else
        x = corrAH(cNew);
    end
    pos=1;
    for c=1:numCondition
        z{c} = x{c}-mean(x{1});
        mData(c,i) = mean(z{c});
        nData(c,i) = length(z{c});
        sData(c,i) = std(z{c})./sqrt(nData(c,i));
        errorbar(hCorrAH,pos,mData(c,i),sData(c,i),'Color',colorName(i,:),'Marker','o');
        if i==2
            if c==1
                text(pos-0.9,mData(c,i)+sData(c,i)+0.005,['N=' num2str(nData(c,i))],'Parent',hCorrAH);
                text(pos-0.9,mData(c,i)-sData(c,i)-0.005,['S=' num2str(nSessions(c))],'Parent',hCorrAH);
            else
                text(pos-0.4,mData(c,i)+sData(c,i)+0.005,num2str(nData(c,i)),'Parent',hCorrAH);
                text(pos-0.2,mData(c,i)-sData(c,i)-0.005,num2str(nSessions(c)),'Parent',hCorrAH);
            end
        end
        hold(hCorrAH,'on')
        
        if mod(c,grp)==0
            pos = pos+2;
        else
            pos=pos+1;
        end
    end
    
    for r=1:3
        k=(r-1)*grp;
        p(i) = plot(hCorrAH,xIndex2(k+1:k+grp),mData(k+1:k+grp,i),'Color',colorName(i,:),'LineWidth',1.2);
    end
    
    
end
xTickRange = [0 numCondition+3];
yTickRange = [-0.08 0.08];
yTick = -0.08:0.02:0.08;
axis(hCorrAH,[xTickRange yTickRange])
set(hCorrAH,'XTick',xIndex2,'XTickLabel',typeList2,'YTick',yTick,'YTicklabel',yTick,'XTickLabelRotation',45,'TickDir','out','box','off')
title(hCorrAH,'Inter-Hemisphere')
xlabel(hCorrAH,'Attention Condition')
ylabel(hCorrAH,'Change in correlation')
legend(hCorrAH,p,'Firing Rate','Gamma Power','Alpha Power','location','northeast')

annotation('textbox',[0.01 0.95 0.02 0.05],'String','\bf A','FontSize',12,'EdgeColor','none')
annotation('textbox',[0.5 0.95 0.02 0.05],'String','\bf B','FontSize',12,'EdgeColor','none')
