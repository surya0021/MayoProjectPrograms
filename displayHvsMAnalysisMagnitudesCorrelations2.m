% This program is used in displayFigure3,displayFigure3v2,displayFigure4,displayFigure4v2
% The difference between this program and displayHvsMAnalysisMagnitudeCorrelations is that here the z-scoring is done relative to the attend out condition of the respective arraySide. 

% dataTypeNum - 1: spikes, 2: gamma, 3: alpha
% processType - 1: raw, 2: z-score

function [dataTMP,mData,mCombinedData,rData,rCombinedData,typeList1,typeList2,nSessions] = displayHvsMAnalysisMagnitudesCorrelations2(dataTypeNum,processType,trialCutoff)

if ~exist('processType','var') || isempty(processType);         processType=2;                  end
if ~exist('trialCutoff','var');         trialCutoff=15;                 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oriChangeList = [2 3]; tpStr = '_TargetOnset'; timePeriod = [-0.5 0]; populationType = 'Stimulated';
tapers = [2 3]; alphaRangeHz = [8 12]; gammaRangeHz = [42 78];

% folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
folderSourceString = 'E:\Mayo';

folderName = fullfile(folderSourceString,'Data','savedDataSummary');
fileNameStr = 'dataOri_';
for i=1:length(oriChangeList)
    fileNameStr = cat(2,fileNameStr,num2str(oriChangeList(i)));
end

fileNameStr = cat(2,fileNameStr,tpStr,num2str(timePeriod(1)),'_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));

fileName = fullfile(folderName,[fileNameStr '.mat']);
load(fileName);

fileNameElectrodes = fullfile(folderName,['electrodeArrayList' populationType '.mat']);
load(fileNameElectrodes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

typeList1 = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'H0N'} {'H1N'} {'M0N'} {'M1N'}];
typeList2 = [{'AIVH'} {'AOVH'} {'AIIH'} {'AOIH'} {'AIVM'} {'AOVM'} {'AIIM'} {'AOIM'} {'TIH'} {'TOH'} {'TIM'} {'TOM'}];

magnitudeTitleStr = [{'Side0'} {'Side1'} {'Combined'}];
correlationsTitleStr = [{'Side0'} {'Side1'} {'Inter-Hemispheric'} {'Combined'}];

colorNameList = 'rbg';
colorName = colorNameList(dataTypeNum);
if dataTypeNum==1
    dataTMP = transformData(spikeData,processType,dataTypeNum,timePeriod,electrodeArrayList);
elseif dataTypeNum==2
    dataTMP = transformData(gammaData,processType,dataTypeNum,timePeriod,electrodeArrayList);
elseif dataTypeNum==3
    dataTMP = transformData(alphaData,processType,dataTypeNum,timePeriod,electrodeArrayList);
end

[numSessions,numConditions] = size(dataTMP{1});
mData = cell(numConditions,2);
rData = cell(numConditions,3);
nTrials = zeros(numSessions,numConditions);
nSessions = zeros(1,numConditions);

for c=1:numConditions
    for s=1:numSessions
            tmp1 = dataTMP{1}{s,c}; %#ok<*USENS>
            tmp2 = dataTMP{2}{s,c};
        
        nTrials(s,c) = size(tmp1,2);
        
        if nTrials(s,c)>trialCutoff
            nSessions(c) = nSessions(c)+1;
            
            % Magnitude
            
            mData{c,1} = cat(1,mData{c,1},mean(tmp1,2));
            mData{c,2} = cat(1,mData{c,2},mean(tmp2,2));
            
            % Correlations
            allCorrs = corrcoef([tmp1' tmp2']);
            n1=size(tmp1,1); n2=size(tmp2,1);
            corr1 = allCorrs(1:n1,1:n1);
            corr2 = allCorrs(n1+(1:n2),n1+(1:n2));
            corr3 = allCorrs(n1+(1:n2),1:n1);
            
            corr1All =corr1(abs(triu(corr1))>0 & abs(triu(corr1))<1);
            rData{c,1} = cat(1,rData{c,1},corr1All);
            corr2All =corr2(abs(triu(corr2))>0 & abs(triu(corr2))<1);
            rData{c,2} = cat(1,rData{c,2},corr2All);
            rData{c,3} = cat(1,rData{c,3},corr3(:));
        else
            disp(['discarded: ' typeList1{c} ', session: ' num2str(s) ', numTrials: ' num2str(nTrials(s,c))]);
        end
    end
end

% Plot Magnitudes
for side=1:2
    h = subplot(3,2,2*(side-1)+1);
    
    mDataTMP = zeros(1,numConditions);
    for c=1:numConditions
        x = mData{c,side};
        mtmp = mean(x);
        ntmp = length(x);
        stmp = std(x)/sqrt(ntmp);
        errorbar(c,mtmp,stmp,[colorName 'o']); hold on;
        text(c-0.2,mtmp+stmp+0.005,num2str(ntmp));
        text(c-0.2,mtmp-stmp-0.005,num2str(nSessions(c)));
        mDataTMP(c) = mtmp;
    end
    plot(mDataTMP,colorName);
    set(h,'XTick',1:12,'XTickLabel',typeList1);
    title(['Magnitude, ' correlationsTitleStr{side}]);
end

h = subplot(3,2,5); % Combine sides
mCombinedData = combineData(mData);

mDataTMP=zeros(1,numConditions);
for c=1:numConditions
    mtmp = mean(mCombinedData{c});
    ntmp = length(mCombinedData{c});
    stmp = std(mCombinedData{c})/sqrt(ntmp);
    
    errorbar(c,mtmp,stmp,[colorName 'o']); hold on;
    text(c-0.2,mtmp+stmp+0.005,num2str(ntmp));
    mDataTMP(c) = mtmp;
end
plot(mDataTMP,colorName);
set(h,'XTick',1:12,'XTickLabel',typeList2);
title(['Magnitude, ' magnitudeTitleStr{3}]);

% Plot Correlations
for side=1:3
    h = subplot(4,2,2*side);
    
    mDataTMP = zeros(1,numConditions);
    for c=1:numConditions
        x = rData{c,side};
        mtmp = mean(x);
        ntmp = length(x);
        stmp = std(x)/sqrt(ntmp);
        errorbar(c,mtmp,stmp,[colorName 'o']); hold on;
        text(c-0.2,mtmp+stmp+0.005,num2str(ntmp));
        text(c-0.2,mtmp-stmp-0.005,num2str(nSessions(c)));
        mDataTMP(c) = mtmp;
    end
    plot(mDataTMP,colorName);
    set(h,'XTick',1:12,'XTickLabel',typeList1);
    title(['Correlations, ' correlationsTitleStr{side}]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = subplot(4,2,8); % Combine sides
rCombinedData = combineData(rData);

mDataTMP=zeros(1,numConditions);
for c=1:numConditions
    mtmp = mean(rCombinedData{c});
    ntmp = length(rCombinedData{c});
    stmp = std(rCombinedData{c})/sqrt(ntmp);
    
    errorbar(c,mtmp,stmp,[colorName 'o']); hold on;
    text(c-0.2,mtmp+stmp+0.005,num2str(ntmp));
    mDataTMP(c) = mtmp;
end
plot(mDataTMP,colorName);
set(h,'XTick',1:12,'XTickLabel',typeList2);
title(['Correlations, ' correlationsTitleStr{4}]);

end

function newData = transformData(data,processType,dataTypeNum,timePeriod,electrodeArrayList)

[numSessions,numConditions,~] = size(data);

if processType==1
    newData{1} = cell(numSessions,numConditions); 
    newData{2} = cell(numSessions,numConditions);
    for s=1:numSessions
        for c=1:numConditions
            for k=1:2
                eList=electrodeArrayList{s}{k};
                if dataTypeNum==1
                    newData{k}{s,c} = cell2mat(squeeze(data(s,c,eList)))./diff(timePeriod);
                else
                    newData{k}{s,c} = cell2mat(squeeze(data(s,c,eList)));
                end
            end
        end
    end
elseif processType==2 % z-score
           
    newData{1} = cell(numSessions,numConditions); 
    newData{2} = cell(numSessions,numConditions);
    for s=1:numSessions
        for k=1:2
            eList=electrodeArrayList{s}{k};
            for e=1:length(eList)
                if dataTypeNum==1
                    x = squeeze(data(s,:,eList(e)));
                    x=cellfun(@(y) y/diff(timePeriod),x,'un',0);
                else
                    x = squeeze(data(s,:,eList(e)));
                end
                allX=[];
                if k==1 % for Right Array
                    c = 2; %H1V
                elseif k==2 % for Left Array
                    c = 1;  %H0V
                end
                allX = cat(2,allX,x{c});
                
                mX = mean(allX);
                stdX = std(allX);
                
                for c=1:numConditions
                    newData{k}{s,c}(e,:) = (x{c}-mX)/stdX;
                end
            end
        end
    end
end
end
function newData = combineData(data)
newData{1} = [data{1,1};data{2,2}]; %AIVH
newData{2} = [data{2,1};data{1,2}]; %AOVH
newData{3} = [data{3,1};data{4,2}]; %AIIH
newData{4} = [data{4,1};data{3,2}]; %AOIH

newData{5} = [data{5,1};data{6,2}]; %AIVM
newData{6} = [data{6,1};data{5,2}]; %AOVM
newData{7} = [data{7,1};data{8,2}]; %AIIM
newData{8} = [data{8,1};data{7,2}]; %AOIM

newData{9} = [data{9,1};data{10,2}]; %TIH
newData{10} = [data{10,1};data{9,2}]; %TOH
newData{11} = [data{11,1};data{12,2}]; %TIM
newData{12} = [data{12,1};data{11,2}]; %TOM
end