function displayAUCMayo(fileNameString,folderSourceString)
if ~exist('folderSourceString','var');  folderSourceString='E:/Mayo';   end

folderName=fullfile(folderSourceString,'Data','savedDataForLDA');
close all;
% Load Data

load(fullfile(folderName,[fileNameString 'All' 'TargetOnset_500ms_Tapers_3']));
class{1}='H0V';   class{2}='H1V';

for i=1:2 % 1=H0V 2=H1V
    FRData{i}=firingRate{2}{i};
    alphaLogPower{i}=log10(alphaData{2}{i});
    gammaLogPower{i}=log10(gammaData{2}{i});
    mFRData{i}=mean(FRData{i},1);
    mAlphaLogPower{i}=mean(alphaLogPower{i},1);
    mGammaLogPower{i}=mean(gammaLogPower{i},1);
    label{i}=repmat({class{i}},size(FRData{i},2),1);
    
end
FRDataAll=[FRData{1} FRData{2}];
mFRDataAll=[mFRData{1} mFRData{2}];

alphaLogPowerAll=[alphaLogPower{1} alphaLogPower{2}];
mAlphaLogPowerAll=[mAlphaLogPower{1} mAlphaLogPower{2}];

gammaLogPowerAll=[gammaLogPower{1} gammaLogPower{2}];
mGammaLogPowerAll=[mGammaLogPower{1} mGammaLogPower{2}];

labelAll=[label{1} ; label{2}]; 

numElec=size(FRData{1},1);

%% Get Area under ROC Curve (AUC)

for k=1:numElec
    aucFR(k)=getAUC(FRDataAll(k,:),labelAll,'H0V');
    aucAlpha(k)=getAUC(alphaLogPowerAll(k,:),labelAll,'H0V');
    aucGamma(k)=getAUC(gammaLogPowerAll(k,:),labelAll,'H0V');
end
meanAUCFR=mean(aucFR);
meanAUCAlpha=mean(aucAlpha);
meanAUCGamma=mean(aucGamma);

aucMeanFR=getAUC(mFRDataAll,labelAll,'H0V');
aucMeanAlpha=getAUC(mAlphaLogPowerAll,labelAll,'H0V');
aucMeanGamma=getAUC(mGammaLogPowerAll,labelAll,'H0V');

%% Plot distribution


numRows=round(sqrt(numElec));
numCol=round(numElec/numRows);
figure(1);
hDist1=getPlotHandles(numRows,numCol,[0.05 0.05 0.9 0.9],0.01,0.03);
figure(2);
hDist2=getPlotHandles(numRows,numCol,[0.05 0.05 0.9 0.9],0.01,0.03);
figure(3);
hDist3=getPlotHandles(numRows,numCol,[0.05 0.05 0.9 0.9],0.01,0.03);
figure(4);
hDist4=getPlotHandles(1,3,[0.1 0.3 0.8 0.5],0.1,0.1);


% Uses Matlab function histfit.m. But does not give option to choose
% binWidth.
% for k=1:numElec
%     axes(hDist(k))
%     h1=histfit(FRData{1}(k,:),10);
%     h1(1).FaceColor=[0 0 1]; h1(2).Color=[0 0 1]; h1(1).FaceAlpha=0.5;
%     hold on
%     h2=histfit(FRData{2}(k,:),10);
%     h2(1).FaceColor=[1 0 0]; h2(2).Color=[1 0 0]; h2(1).FaceAlpha=0.2;
%     title(['elec ' num2str(k)],'FontSize',8)
%   xlim([0 70])
% end

% Firing Rate
for m=1:numElec
    makeHist(hDist1(m),FRData{1}(m,:),[0 0 1],0.5,5);
    hold on
    makeHist(hDist1(m),FRData{2}(m,:),[1 0 0],0.2,5);
    xlim(hDist1(m),[0 100]);
    text(0.5,0.6,['auc=' num2str(aucFR(m))],'Units','normalized','FontSize',8,'Parent',hDist1(m));
    title(hDist1(m),['elec ' num2str(m)],'FontSize',8);
end
axes(hDist1)
set(gcf,'defaultaxesfontsize',8);
suptitle(['Firing Rate: mAUC=' num2str(meanAUCFR)]);
% Alpha Power
for m=1:numElec
    makeHist(hDist2(m),alphaLogPower{1}(m,:),[0 0 1],0.5,0.1);
    hold on
    makeHist(hDist2(m),alphaLogPower{2}(m,:),[1 0 0],0.2,0.1);
    text(0.5,0.6,['auc=' num2str(aucAlpha(m))],'Units','normalized','FontSize',8,'Parent',hDist2(m));
    title(hDist2(m),['elec ' num2str(m)],'FontSize',8);
end
axes(hDist2)
set(gcf,'defaultaxesfontsize',8);
suptitle(['Alpha Power: mAUC=' num2str(meanAUCAlpha)]);
%Gamma Power
for m=1:numElec
    makeHist(hDist3(m),gammaLogPower{1}(m,:),[0 0 1],0.5,0.1);
    hold on
    makeHist(hDist3(m),gammaLogPower{2}(m,:),[1 0 0],0.2,0.1);
    text(0.5,0.6,['auc=' num2str(aucGamma(m))],'Units','normalized','FontSize',8,'Parent',hDist3(m));
    title(hDist3(m),['elec ' num2str(m)],'FontSize',8);
end
axes(hDist3)
set(gcf,'defaultaxesfontsize',8);
suptitle(['Gamma Power: mAUC=' num2str(meanAUCGamma)]);

% Mean Firing Rate
makeHist(hDist4(1),mFRData{1},[0 0 1],0.5,5);
hold on
makeHist(hDist4(1),mFRData{2},[1 0 0],0.2,5);
title(hDist4(1),['Firing Rate: AUC=' num2str(aucMeanFR)]);
% Mean Alpha Power
makeHist(hDist4(2),mAlphaLogPower{1},[0 0 1],0.5,0.1);
hold on
makeHist(hDist4(2),mAlphaLogPower{2},[1 0 0],0.2,0.1);
title(hDist4(2),['Alpha Power: AUC=' num2str(aucMeanAlpha)]);
%Mean Gamma Power
makeHist(hDist4(3),mGammaLogPower{1},[0 0 1],0.5,0.1);
hold on
makeHist(hDist4(3),mGammaLogPower{2},[1 0 0],0.2,0.1);
title(hDist4(3),['Gamma Power: AUC=' num2str(aucMeanGamma)]);


end


% makeHist function plots the histogram and fits the gaussian
% distribution to the histogram. Unlike histfit function of Matlab,one can
% choose binWidth.(Built in line with Matlab's histfit function)
function makeHist(h,data,Color,FA,BW) % FA=FaceAlpha %BW=BinWidth
n=numel(data);
[counts,edges]=histcounts(data,'BinWidth',BW);
centers = edges(1:end-1)+diff(edges)/2;
bar(h,centers,counts,1,'FaceColor',Color,'FaceAlpha',FA);

pd=fitdist(data','normal');
% Find range for plotting 
q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
x = linspace(q(1),q(2));

% Normalize the density to match the total area of the histogram
width = edges(2)-edges(1); % Finds the width of each bin
area = n * width;
y = area * pdf(pd,x);

set(h,'NextPlot','add')
plot(h,x,y,'color',Color,'LineWidth',2)

end


function auc=getAUC(data,label,posClass)

posClassData=data(find(strcmp(label,posClass)));
negClassData=data(find(~strcmp(label,posClass)));
criterion=linspace(min(data),max(data),50);

for i=1:length(criterion)
    TP(i)=numel(find(posClassData >= criterion(i))); %True Positive
    TN(i)=numel(find(negClassData < criterion(i))); %True Negative
    FP(i)=numel(find(negClassData >= criterion(i))); %False Positive
    FN(i)=numel(find(posClassData < criterion(i))); %False negative
end
TPR= TP./(TP+FN); % True Positive Rate or Recall
FPR= FP./(FP+TN); % False Postive Rate
auc=abs(trapz(FPR,TPR));
auc=abs(auc-0.5);
end
    
