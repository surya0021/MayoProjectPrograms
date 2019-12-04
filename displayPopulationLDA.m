% This function saves the LDA or Regularized LDA projections of 12 cue conditions for different population size / features and Gamma
% Parameter

% dataType: 1:Firing Rate; 2:Gamma Power; 3: Alpha Power; 4:Firing Rate + Gamma Power; 5:Firing Rate + Gamma Power + Alpha Power   
% regFlag: 0: LDA ;      1: Regularized LDA ; 
function displayPopulationLDA(dataType,regFlag) 

if ~exist('regFlag','var');     regFlag=0;          end
oriChangeList = [2 3]; tpStr = '_TargetOnset'; timePeriod = [-0.5 0]; populationType = 'Stimulated';
tapers = [2 3];
alphaRangeHz = [8 12]; gammaRangeHz = [42 78];  trialCutoff = 15;

folderSourceString = 'E:\Mayo';
folderNameData = fullfile(folderSourceString,'Data','savedDataSummary');
folderNameSave = fullfile(folderSourceString,'Data','populationSizeLDA');
fileNameStr = 'dataOri_';

for i=1:length(oriChangeList)
    fileNameStr = cat(2,fileNameStr,num2str(oriChangeList(i)));
end

fileNameStr = cat(2,fileNameStr,tpStr,num2str(timePeriod(1)),'_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));

fileNameData = fullfile(folderNameData,[fileNameStr '.mat']);
load(fileNameData);

fileNameSaveElectrodes = fullfile(folderNameData,['electrodeArrayList' populationType '.mat']);
load(fileNameSaveElectrodes);
eList = cellfun(@(x) cat(1,x{1},x{2}),electrodeArrayList,'un',0);
numFolds = 5;
if regFlag==0
    ldaStr='LDA';
else
    ldaStr='RLDA';
end
if dataType==1
    data = cellfun(@(x) x./diff(timePeriod),spikeData,'un',0); % converting spike count data into firing rate
    dataStr = 'FR';
elseif dataType==2
    data = gammaData;
    dataStr = 'Gamma';
elseif dataType==3
    data = alphaData;
    dataStr = 'Alpha';
elseif dataType==4
    data{1} = cellfun(@(x) x./diff(timePeriod),spikeData,'un',0); % converting spike count data into firing rate
    data{2} = gammaData;
    dataStr = 'FR + Gamma';
elseif dataType==5
    data{1} = cellfun(@(x) x./diff(timePeriod),spikeData,'un',0); % converting spike count data into firing rate
    data{2} = gammaData;
    data{3} = alphaData;
    dataStr = 'FRGammaAlpha';
    
end

disp(ldaStr)
if dataType>3
    [nSessions,nConditions,~] = size(data{1});
else
    [nSessions,nConditions,~] = size(data);
end
if regFlag==0
    [projection,Gamma] = getProjection(data,eList,numFolds,trialCutoff,regFlag,dataType); %#ok<ASGLU>
else
    [projection,Gamma,gammaSampled,errors] = getProjection(data,eList,numFolds,trialCutoff,regFlag,dataType); %#ok<ASGLU>
end
meanProj = cellfun(@(y) cell2mat(cellfun(@mean,cellfun(@(x) mean(x,2),y,'un',0),'un',0)),projection,'un',0); % mean across electrode set iterations and cv folds
maxPoolSize = max(cellfun(@length,meanProj(:,1)));
mProj=cell(12,maxPoolSize);

for poolSize=1:maxPoolSize
    for c=1:nConditions
        for s=1:nSessions
            if length(meanProj{s,c})>=poolSize
                if ~isnan(meanProj{s,c}(poolSize))
                    mProj{c,poolSize} = cat(1,mProj{c,poolSize},meanProj{s,c}(poolSize));
                else
                    continue
                end
            else
                continue
            end
        end
    end
end
fileNameStr2 = cat(2,[dataStr,ldaStr],'_','Folds',num2str(numFolds),tpStr,num2str(timePeriod(1)), '_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));

fileNameSave=fullfile(folderNameSave,[fileNameStr2,'.mat']);

if regFlag==0
    save(fileNameSave,'mProj','Gamma')
else
    save(fileNameSave,'mProj','Gamma','gammaSampled','errors')
end

end


function [proj,paramGamma,gammaSampled,errors] = getProjection(data,eList,numFolds,trialCutoff,regFlag,dataType)
if dataType>3
    [nSessions,nConditions,~] = size(data{1});
else
    [nSessions,nConditions,~] = size(data);
end
proj = cell(nSessions,nConditions);
for s=1:nSessions
    disp(['Working on session ' num2str(s) ' of ' num2str(nSessions)])
    nSet = length(eList{s});
    for c=1:nConditions
        proj{s,c} = cell(1,nSet);
    end
    eSet = getElectrodeSet(eList{s},nSet);
    if dataType>3
        label1 = repmat({'H0V'},size(data{1}{s,1,1},2),1);
        label2 = repmat({'H1V'},size(data{1}{s,2,1},2),1);
    else
        label1 = repmat({'H0V'},size(data{s,1,1},2),1);
        label2 = repmat({'H1V'},size(data{s,2,1},2),1);
    end
    labelAll = cat(1,label1,label2);
    class = unique(labelAll);
    cvp = cvpartition(labelAll,'KFold',numFolds);

    for poolSize=1:length(eList{s})
        disp(['PoolSize ' num2str(poolSize) ' of ' num2str(length(eList{s}))])
        for nIter= 1:size(eSet{poolSize},1)
            if dataType>3
                data1 = [];   data2 = [];   
                for p=1:length(data)
                    data1 = cat(1,data1,cell2mat(squeeze(data{p}(s,1,eSet{poolSize}(nIter,:))))); %H0V
                    data2 = cat(1,data2,cell2mat(squeeze(data{p}(s,2,eSet{poolSize}(nIter,:))))); %H1V
                end
            else
                data1 = cell2mat(squeeze(data(s,1,eSet{poolSize}(nIter,:)))); %H0V
                data2 = cell2mat(squeeze(data(s,2,eSet{poolSize}(nIter,:)))); %H1V
            end
            dataAll = cat(2,data1,data2);
            cvPos = ones(1,12);
            for cv=1:numFolds
                
                traininds=cvp.training(cv);
                testinds=cvp.test(cv);
                Mdl = fitcdiscr(dataAll(:,traininds)',labelAll(traininds));
                if regFlag==1
                    Mdl = fitcdiscr(dataAll',labelAll);
                    [err,gamma,~,~]=cvshrink(Mdl,'NumGamma',30);
                    [~,minidx]=min(err);
                    newGamma=gamma(minidx);
                    Mdl.Gamma=newGamma;
                    errors{s}{poolSize}(:,nIter,cv) = err;
                    gammaSampled{s}{poolSize}(:,nIter,cv) = gamma;   
                end
                Mu1=Mdl.Mu(1,:)';
                Mu2=Mdl.Mu(2,:)';
                SW=Mdl.Sigma;
                paramGamma{s}{poolSize}(nIter,cv) = Mdl.Gamma;
                if isrow(SW) && size(dataAll,1)>1    % For gamma=1 case SW is a row matrix (diag(s.^2))
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
                
                for c=1:nConditions
                    if c==1 || c==2
                        dataTMP = dataAll(:,testinds);
                        dataTMP = dataTMP(:,strcmp(labelAll(testinds),class{c}));
                    else
                        if dataType>3
                            dataTMP=[];
                            for p=1:length(data)
                                dataTMP = cat(1,dataTMP,cell2mat(squeeze(data{p}(s,c,eSet{poolSize}(nIter,:)))));
                            end
                        else
                            dataTMP = cell2mat(squeeze(data(s,c,eSet{poolSize}(nIter,:))));
                        end
                    end
                    newData{c} = dataTMP'*weightVector;
                end
                mX = mean(newData{2});
                stdX = std(newData{2});
                for c=1:nConditions
                    z{c} = (newData{c}-mX)/stdX;
                    if length(z{c})>trialCutoff
                        proj{s,c}{poolSize}(nIter,cvPos(c))=mean(z{c});
                        cvPos(c)=cvPos(c)+1;
                    else
                        continue
                    end
                end
%                 proj{s}{poolSize}{nIter,cv}=z;
                clear z 
            end            
        end
    end
end
end

function eSet = getElectrodeSet(eList,nSet)
for i=1:length(eList)
    if i==1
        eSet{i} = eList; %#ok<*AGROW>
    elseif i==length(eList)
        eSet{i} = eList';
    else
        n=1;
        while n<=nSet
            eSet{i}(n,:) = sort(randsample(eList,i));
            eSet{i} = unique(eSet{i},'rows','stable');
            n=size(eSet{i},1)+1;
        end
    end
end
end