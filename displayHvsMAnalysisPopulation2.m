% This Program is used in displayFigure5, displayFigure5v2 to save the
% output of this function. 
% The difference between this function and
% displayHvsMPopulation is that here pooled variance is used for uncorrelated LDA case 

% dataTypeNum - 1: spikes, 2: gamma, 3: alpha
% transformType - 1: simple averaging across sides, 2: LDA-uncorrelated,
% 3: LDA-covariance

function [mData,typeList1] = displayHvsMAnalysisPopulation2(dataTypeNum,transformType,trialCutoff,normalizeFlag)

if ~exist('transformType','var');       transformType=2;                end
if ~exist('trialCutoff','var');         trialCutoff=15;                 end
if ~exist('normalizeFlag','var');       normalizeFlag=1;                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oriChangeList = [2 3]; tpStr = '_TargetOnset'; timePeriod = [-0.5 0]; populationType = 'Stimulated';
tapers = [2 3]; alphaRangeHz = [8 12]; gammaRangeHz = [42 78];

% folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
folderSourceString = 'E:\Mayo';
folderNameSave = fullfile(folderSourceString,'Data','savedDataSummary');
fileNameStr = 'dataOri_';
for i=1:length(oriChangeList)
    fileNameStr = cat(2,fileNameStr,num2str(oriChangeList(i)));
end

fileNameStr = cat(2,fileNameStr,tpStr,num2str(timePeriod(1)),'_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));

fileNameSave = fullfile(folderNameSave,[fileNameStr '.mat']);
load(fileNameSave);

fileNameSaveElectrodes = fullfile(folderNameSave,['electrodeArrayList' populationType '.mat']);
load(fileNameSaveElectrodes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

typeList1 = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'H0N'} {'H1N'} {'M0N'} {'M1N'}];

colorNameList = 'rbg';
colorName = colorNameList(dataTypeNum);
if dataTypeNum==1
    dataTMP = transformData(spikeData,transformType,dataTypeNum,timePeriod,electrodeArrayList);
elseif dataTypeNum==2
    dataTMP = transformData(gammaData,transformType,dataTypeNum,timePeriod,electrodeArrayList);
elseif dataTypeNum==3
    dataTMP = transformData(alphaData,transformType,dataTypeNum,timePeriod,electrodeArrayList);
end

if normalizeFlag
    dataTMP = normalizeData(dataTMP);
end

[numSessions,numConditions] = size(dataTMP);
mData = cell(1,numConditions);
nTrials = zeros(numSessions,numConditions);

for c=1:numConditions
    pos=1;
    for s=1:numSessions
        tmp = dataTMP{s,c};
        nTrials(s,c) = length(tmp);
        if nTrials(s,c)>trialCutoff
            mData{c}(pos) = mean(tmp);
            pos = pos+1;
        else
            continue
        end
    end
end

% Plot Data
mDataTMP = zeros(1,numConditions);
for c=1:numConditions
    
    x = mData{c};
    mtmp = mean(x);
    ntmp = length(x);
    stmp = std(x)/sqrt(ntmp);
    errorbar(c,mtmp,stmp,[colorName 'o']); hold on;
    text(c-0.2,mtmp+stmp+0.005,num2str(ntmp));
    mDataTMP(c) = mtmp;
end
plot(mDataTMP,colorName);
h=gca;
set(h,'XTick',1:12,'XTickLabel',typeList1);

end

function newData = transformData(data,transformType,dataTypeNum,timePeriod,electrodeArrayList)

[numSessions,numConditions,~] = size(data);
newData=cell(numSessions,numConditions);
    
for s=1:numSessions
    eList = [electrodeArrayList{s}{1};electrodeArrayList{s}{2}]';
    if dataTypeNum==1
        data1 = cell2mat(squeeze(data(s,1,eList)))./diff(timePeriod); % H0V
        data2 = cell2mat(squeeze(data(s,2,eList)))./diff(timePeriod); % H1V
    else
        data1 = cell2mat(squeeze(data(s,1,eList)));
        data2 = cell2mat(squeeze(data(s,2,eList)));
    end
    % Get the weightVector
    if transformType==1 % Simple Averaging of sides
        n1=length(electrodeArrayList{s}{1});
        n2=length(electrodeArrayList{s}{2});
        weightVector = [ones(n1,1)/n1; -ones(n2,1)/n2];
        
    elseif transformType==2 % LDA for uncorrelated case
        m1 = size(data1,2);
        m2= size(data2,2);
        meanDiff = mean(data1,2) - mean(data2,2);
        var1 = var(data1,[],2);
        var2 = var(data2,[],2);
        pooledVar = ((m1-1)*var1 + (m2-1)*var2)/(m1+m2-2); %Pooled Variance
%       stdData = sqrt((var(data1,[],2)+var(data2,[],2))/2); % Average std
%       stdData = sqrt(pooledVar); % Pooled standard deviation
        weightVector = meanDiff./pooledVar;
        
    elseif transformType==3 % LDA
        label1 = repmat({'H0V'},size(data1,2),1);
        label2 = repmat({'H1V'},size(data2,2),1);
        labelAll = cat(1,label1,label2);
        dataAll = cat(2,data1,data2);
        
        % LDA using fitcdiscr function of MATLAB
  
        Mdl = fitcdiscr(dataAll',labelAll);
        Mu1 = Mdl.Mu(1,:)';
        Mu2 = Mdl.Mu(2,:)';
        SW = Mdl.Sigma;
        if isrow(SW) && size(data,1)>1    % For gamma=1 case SW is a row matrix (diag(s.^2))
            SW=bsxfun(@times,SW,eye(length(SW))); % converts into a diagonal matrix
        end
        
        SB = (Mu1-Mu2)*(Mu1-Mu2)'; % between class scatter
        if det(SW)==0
            error('Within scatter matrix is non-invertible')
        end
        if rank(SB)~=1
            error('The rank of SB is more than 1')
        end
        weightVector = SW\(Mu1-Mu2);
    end
     
    for c=1:numConditions
        if dataTypeNum==1
            dataTMP = cell2mat(squeeze(data(s,c,eList)))./diff(timePeriod);
        else
            dataTMP = cell2mat(squeeze(data(s,c,eList)));
        end
        newData{s,c} = dataTMP' * weightVector;
    end
end
end
function newData = normalizeData(data)

[numSessions,numConditions] = size(data);
newData=cell(numSessions,numConditions);

for s=1:numSessions 
    mX = mean(data{s,2});   % Normalized w.r.t H1V condition
%   stdX = sqrt((var(data{s,1}) + var(data{s,2}))/2);
    stdX= std(data{s,2});

        
    for c=1:numConditions
        newData{s,c} = (data{s,c}-mX)/stdX;
    end
end
end