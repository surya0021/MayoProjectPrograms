SessionID = 'all'; 
neuronType = 'All';
populationType = 'Stimulated';
OriChange = [2 3];  
tapers = [2 3];
folderSourceString = 'E:/Mayo';
folderSave = fullfile(folderSourceString,'Data','FigureData');
fileToSave  = fullfile(folderSave,['Figure5Dataset_' SessionID neuronType populationType 'Ori' num2str(OriChange(1)) num2str(OriChange(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '.mat']);
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
        Figure5Dataset(dataTypeNum) = struct('dataType',dataType{dataTypeNum},'meanDiff', mData(1),...
            'uncorrelatedLDAData',mData(2),'LDA',mData(3),'attCueList',{typeList1}); %#ok<SAGROW>
        
    end
    save(fileToSave,'Figure5Dataset')
end
colorName = gray(4);
attCueList = Figure5Dataset(1).attCueList;
typeList1 = [{'H1V'} {'M1V'} {'M0V'} {'H0V'} {'H1N'} {'M1N'} {'M0N'} {'H0N'} {'M1I'} {'H1I'} {'H0I'} {'M0I'}];
cNew = [2 6 5 1 10 12 11 9 8 4 3 7];
xIndex = [1:4 6:9 11:14];

%Figure 5

set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]) % Figure 4
hPlot = getPlotHandles(3,1,[0.07 0.1 0.9 0.85],0.03,0.02,0);  linkaxes(hPlot);
for i=1:3
    meanDiff = Figure5Dataset(i).meanDiff(cNew);
    unCorrLDA = Figure5Dataset(i).uncorrelatedLDAData(cNew);
    LDA = Figure5Dataset(i).LDA(cNew);
    mDataMean = cell2mat(cellfun(@mean,meanDiff,'un',0));
    nDataMean= cell2mat(cellfun(@length,meanDiff,'Un',0));
    mDataUncorr = cell2mat(cellfun(@mean,unCorrLDA,'un',0));
    mDataLDA = cell2mat(cellfun(@mean,LDA,'un',0));
    sDataMean = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),meanDiff,'un',0));
    sDataUncorr = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),unCorrLDA,'un',0));
    sDataLDA = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),LDA,'un',0));
    r=1; pos=1;
    for c = 1:length(mDataMean)
        errorbar(hPlot(1),pos,mDataMean(c),sDataMean(c),'Color',colorName(i,:),'Marker','o'); hold(hPlot(1),'on');
        errorbar(hPlot(2),pos,mDataUncorr(c),sDataUncorr(c),'Color',colorName(i,:),'Marker','o'); hold(hPlot(2),'on')
        errorbar(hPlot(3),pos,mDataLDA(c),sDataLDA(c),'Color',colorName(i,:),'Marker','o'); hold(hPlot(3),'on')
        if i==1
            if c==1
                text(pos-0.3,mDataMean(c)+sDataMean(c)+0.1,['N=' num2str(nDataMean(c))],'Parent',hPlot(1));
            else
                text(pos-0.1,mDataMean(c)+sDataMean(c)+0.2,num2str(nDataMean(c)),'Parent',hPlot(1));
            end
        end
        
        if mod(c,4)==0
           pos = pos+2;
        else
           pos = pos+1; 
        end
    end
    for r=1:3
        k=(r-1)*4;
        p1(i) = plot(hPlot(1),xIndex(k+1:k+4),mDataMean(k+1:k+4),'Color',colorName(i,:),'LineWidth',1.2);   hold(hPlot(1),'on')
        plot(hPlot(2),xIndex(k+1:k+4),mDataUncorr(k+1:k+4),'Color',colorName(i,:),'LineWidth',1.2);     hold(hPlot(2),'on')
        plot(hPlot(3),xIndex(k+1:k+4),mDataLDA(k+1:k+4),'Color',colorName(i,:),'LineWidth',1.2);    hold(hPlot(3),'on')
    end
    line([0 15],[mDataMean(4) mDataMean(4)],'Color',colorName(i,:),'Parent',hPlot(1))
    line([0 15],[mDataUncorr(4) mDataUncorr(4)],'Color',colorName(i,:),'Parent',hPlot(2))
    line([0 15],[mDataLDA(4) mDataLDA(4)],'Color',colorName(i,:),'Parent',hPlot(3))
end
line([0 15],[mDataMean(1) mDataMean(1)],'Color',colorName(1,:),'LineStyle','--','Parent',hPlot(1))
line([0 15],[mDataUncorr(1) mDataUncorr(1)],'Color',colorName(1,:),'LineStyle','--','Parent',hPlot(2))
line([0 15],[mDataLDA(1) mDataLDA(1)],'Color',colorName(1,:),'LineStyle','--','Parent',hPlot(3))

%%%%%%%% Axis Configuration %%%%%%%%%%%%%
fontSize = 10;
tickLengthPlot = get(hPlot(1),'TickLength'); 
xRange = [0 15];
yRange = [-0.5 1.5];
yTickRange = -0.5:0.5:1.5;
for i=1:3
    axis(hPlot(i),[xRange yRange])
end
for r=1:3
    k=(r-1)*4;
    set(hPlot(1),'XTick',xIndex,'XTickLabel',[],'XTickLabelRotation',45,'Ytick',yTickRange,'YTicklabel',yTickRange,'fontSize',fontSize,'TickDir','out','box','off');
    set(hPlot(2),'XTick',xIndex,'XTickLabel',[],'XTickLabelRotation',45,'Ytick',yTickRange,'YTicklabel',yTickRange,'fontSize',fontSize,'TickDir','out','box','off');
    set(hPlot(3),'XTick',xIndex,'XTickLabel',typeList1,'XTickLabelRotation',45,'Ytick',yTickRange,'YTicklabel',yTickRange,'fontSize',fontSize,'TickDir','out','box','off');
end
legend(hPlot(1),p1,'Firing Rate','Gamma Power','Alpha Power','location','northwest');
xlabel(hPlot(3),'Attention Condition')
ylabel(hPlot(1),{'\bf Mean Difference';'\rm Z-Scored Projection'})
ylabel(hPlot(2),{'\bf Uncorrelated LDA';'\rm Z-Scored Projection'})
ylabel(hPlot(3),{'\bf LDA';'\rm Z-Scored Projection'})


