SessionID = 'all'; 
neuronType = 'All';
populationType = 'Stimulated';
OriChange = [2 3];  
tapers = [2 3];
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

hCorrWH = getPlotHandles(1,3,[0.05 0.1 0.4 0.85],0.03,0.02,0);  linkaxes(hCorrWH)
hCorrAH = getPlotHandles(1,3,[0.55 0.1 0.4 0.85],0.03,0.02,0);  linkaxes(hCorrAH)

attCueList = Figure3Dataset(1).attCueList;
nSessions = Figure3Dataset(1).nSessions;
attCueListCombined = Figure3Dataset(1).attCueListCombined;
typeList1 = [{'AOVH'} {'AOVM'} {'AIVM'} {'AIVH'} {'TOH'} {'TOM'} {'TIM'} {'TIH'} {'AOIM'} {'AOIH'} {'AIIH'} {'AIIM'}];
typeList2 = [{'H1V'} {'M1V'} {'M0V'} {'H0V'} {'H1N'} {'M1N'} {'M0N'} {'H0N'} {'M1I'} {'H1I'} {'H0I'} {'M0I'}];
cNew = [2 6 5 1 10 12 11 9 8 4 3 7];
nSessions = nSessions(cNew);
colorName = 'rbg';

% Intra-Hemispheric Correlation
mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);   
for i=1:3
    x = Figure3Dataset(i).correlation_WH;
    x = x(cNew);
    r=1; pos=1;
    for c=1:12
        z{c} = x{c}-mean(x{1}); %#ok<SAGROW>
        mData(c,i) = mean(z{c});
        nData(c,i) = length(z{c});
        sData(c,i) = std(z{c})./sqrt(nData(c,i));
        errorbar(hCorrWH(r),pos,mData(c,i),sData(c,i),[colorName(i),'o']);
        if i==2
            if c==1
                text(pos-0.9,mData(c,i)+sData(c,i)+0.005,['N=' num2str(nData(c,i))],'Parent',hCorrWH(r));
                text(pos-0.9,mData(c,i)-sData(c,i)-0.005,['S=' num2str(nSessions(c))],'Parent',hCorrWH(r));
            else
                text(pos-0.4,mData(c,i)+sData(c,i)+0.005,num2str(nData(c,i)),'Parent',hCorrWH(r));
                text(pos-0.2,mData(c,i)-sData(c,i)-0.005,num2str(nSessions(c)),'Parent',hCorrWH(r));
            end
        end
        hold(hCorrWH(r),'on')
        pos=pos+1;
        if mod(c,4)==0
            r = r+1; pos = 1;
        end
    end
    for r=1:3
        k=(r-1)*4;
        p(i) = plot(hCorrWH(r),1:4,mData(k+1:k+4,i),colorName(i));
    end
    
end
YTickCWHRange = [-0.1 0.1];
YTickCWH = -0.1:0.02:0.1;
for r=1:3
    k=(r-1)*4;
    axis(hCorrWH(1),[[0 5] YTickCWHRange]);
%     set(hCorrWH(r),'XTick',1:4,'XTickLabel',[],'YTick',YTickCWH,'YTickLabel',YTickCWH,'TickDir','out','box','off');
    set(hCorrWH(r),'XTick',1:4,'XTickLabel',typeList1(k+1:k+4),'YTick',YTickCWH,'YTickLabel',YTickCWH,'XTickLabelRotation',45,'TickDir','out','box','off');
end

title(hCorrWH(2),'Intra-Hemisphere')
xlabel(hCorrWH(2),'Attention Condition')
ylabel(hCorrWH(1),'Change in correlation')
legend(hCorrWH(1),p,'Firing Rate','Gamma Power','Alpha Power','location','northwest')

% Inter-Hemispheric Correlation
mData = zeros(12,3); sData = zeros(12,3);   nData = zeros(12,3);   
for i=1:3
    x = Figure3Dataset(i).correlation_AH;
    x = x(cNew);
    r=1; pos=1;
    for c=1:12
        z{c} = x{c}-mean(x{1});
        mData(c,i) = mean(z{c});
        nData(c,i) = length(z{c});
        sData(c,i) = std(z{c})./sqrt(nData(c,i));
        errorbar(hCorrAH(r),pos,mData(c,i),sData(c,i),[colorName(i),'o']);
        if i==2
            if c==1
                text(pos-0.9,mData(c,i)+sData(c,i)+0.005,['N=' num2str(nData(c,i))],'Parent',hCorrAH(r));
                text(pos-0.9,mData(c,i)-sData(c,i)-0.005,['S=' num2str(nSessions(c))],'Parent',hCorrAH(r));
            else
                text(pos-0.4,mData(c,i)+sData(c,i)+0.005,num2str(nData(c,i)),'Parent',hCorrAH(r));
                text(pos-0.2,mData(c,i)-sData(c,i)-0.005,num2str(nSessions(c)),'Parent',hCorrAH(r));
            end
        end
        hold(hCorrAH(r),'on')
        pos=pos+1;
        if mod(c,4)==0
            r = r+1; pos = 1;
        end
    end
    for r=1:3
        k=(r-1)*4;
        p(i) = plot(hCorrAH(r),1:4,mData(k+1:k+4,i),colorName(i));
    end
    
end
xTickRange = [0 5];
yTickRange = [-0.1 0.1];
yTick = -0.1:0.02:0.1;
axis(hCorrAH(1,r),[xTickRange yTickRange])
for r=1:3
    k=(r-1)*4;
%     set(hCorrAH(r),'XTick',1:4,'XTickLabel',[],'YTick',yTick,'YTicklabel',yTick,'TickDir','out','box','off');
    set(hCorrAH(r),'XTick',1:4,'XTickLabel',typeList2(k+1:k+4),'YTick',yTick,'YTicklabel',yTick,'XTickLabelRotation',45,'TickDir','out','box','off');
end

title(hCorrAH(2),'Inter-Hemisphere')
xlabel(hCorrAH(2),'Attention Condition')
ylabel(hCorrAH(1),'Change in correlation')
legend(hCorrAH(1),p,'Firing Rate','Gamma Power','Alpha Power','location','northwest')



