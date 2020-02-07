SessionID = 'all'; 
neuronType = 'All';
populationType = 'Stimulated';
OriChange = [2 3];  
tapers = [2 3];
folderSourceString = 'E:/Mayo';
folderSave = fullfile(folderSourceString,'Data','FigureData');
fileToSave  = fullfile(folderSave,['Figure4Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '.mat']);
dataType{1} = 'Firing Rate';
dataType{2} = 'Gamma Power';
dataType{3} = 'Alpha Power';

if exist(fileToSave,'file')
    load(fileToSave)
else
    for dataTypeNum = 1:3
        mData = cell(1,3);
        disp(['Working on ' dataType{dataTypeNum}]);
        for transformType=1:3
            [mData{transformType},typeList1] = displayHvsMAnalysisPopulation2(dataTypeNum,transformType);
        end
        Figure4Dataset(dataTypeNum) = struct('dataType',dataType{dataTypeNum},'meanDiff', mData(1),...
            'uncorrelatedLDAData',mData(2),'LDA',mData(3),'attCueList',{typeList1}); %#ok<SAGROW>
        
    end
    save(fileToSave,'Figure4Dataset')
end
colorName = 'rbg';
attCueList = Figure4Dataset(1).attCueList;
typeList1 = [{'H1V'} {'M1V'} {'M0V'} {'H0V'} {'H1N'} {'M1N'} {'M0N'} {'H0N'} {'M1I'} {'H1I'} {'H0I'} {'M0I'}];
cNew = [2 6 5 1 10 12 11 9 8 4 3 7];

%Figure 4

set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 4

for i=1:3
    meanDiff = Figure4Dataset(i).meanDiff(cNew);
    unCorrLDA = Figure4Dataset(i).uncorrelatedLDAData(cNew);
    LDA = Figure4Dataset(i).LDA(cNew);
    hPlot = getPlotHandles(3,3,[0.07 0.1 0.9 0.85],0.03,0.02,0);  linkaxes(hPlot);
    mDataMean = cell2mat(cellfun(@mean,meanDiff,'un',0));
    nDataMean= cell2mat(cellfun(@length,meanDiff,'Un',0));
    mDataUncorr = cell2mat(cellfun(@mean,unCorrLDA,'un',0));
    mDataLDA = cell2mat(cellfun(@mean,LDA,'un',0));
    sDataMean = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),meanDiff,'un',0));
    sDataUncorr = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),unCorrLDA,'un',0));
    sDataLDA = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),LDA,'un',0));
    r=1; pos=1;
    for c = 1:length(mDataMean)
        errorbar(hPlot(i,r),pos,mDataMean(c),sDataMean(c),[colorName(1),'o']); hold(hPlot(i,r),'on');
        errorbar(hPlot(i,r),pos,mDataUncorr(c),sDataUncorr(c),[colorName(2),'o']);
        errorbar(hPlot(i,r),pos,mDataLDA(c),sDataLDA(c),[colorName(3),'o']);
        if i==1
            if c==1
                text(pos-0.3,mDataLDA(c)+sDataLDA(c)+0.1,['N=' num2str(nDataMean(c))],'Parent',hPlot(i,r));
            else
                text(pos-0.1,mDataLDA(c)+sDataLDA(c)+0.07,num2str(nDataMean(c)),'Parent',hPlot(i,r));
            end
        end
        pos=pos+1;
        if mod(c,4)==0
            r = r+1; pos = 1;
        end
    end
    for r=1:3
        k=(r-1)*4;
        p1 = plot(hPlot(i,r),1:4,mDataMean(k+1:k+4),colorName(1)); hold(hPlot(i,r),'on')
        p2 = plot(hPlot(i,r),1:4,mDataUncorr(k+1:k+4),colorName(2));
        p3 = plot(hPlot(i,r),1:4,mDataLDA(k+1:k+4),colorName(3));
    end
end

%%%%%%%% Axis Configuration %%%%%%%%%%%%%
fontSize = 10;
tickLengthPlot = get(hPlot(1),'TickLength'); 
xRange = [0 5];
yRange = [-0.5 1.5];
yTickRange = -0.5:0.5:1.5;
for i=1:3
    axis(hPlot(i),[xRange yRange])
end
for r=1:3
    k=(r-1)*4;
    set(hPlot(1,r),'XTick',1:4,'XTickLabel',[],'XTickLabelRotation',45,'Ytick',yTickRange,'YTicklabel',yTickRange,'fontSize',fontSize,'TickDir','out','box','off');
    set(hPlot(2,r),'XTick',1:4,'XTickLabel',[],'XTickLabelRotation',45,'Ytick',yTickRange,'YTicklabel',yTickRange,'fontSize',fontSize,'TickDir','out','box','off');
    set(hPlot(3,r),'XTick',1:4,'XTickLabel',typeList1(k+1:k+4),'XTickLabelRotation',45,'Ytick',yTickRange,'YTicklabel',yTickRange,'fontSize',fontSize,'TickDir','out','box','off');
end
legend(hPlot(1),[p1 p2 p3],'Mean Difference','UncorrelatedLDA','LDA','location','northwest');
xlabel(hPlot(3),'Attention Condition')
ylabel(hPlot(1,1),{'\bf Firing Rate';'\rm Z-Scored Projection'})
ylabel(hPlot(2,1),{'\bf Gamma Power';'\rm Z-Scored Projection'})
ylabel(hPlot(3,1),{'\bf Alpha Power';'\rm Z-Scored Projection'})


