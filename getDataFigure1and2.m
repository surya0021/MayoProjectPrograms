% This program generates data to make figure 1(PSD and PSTH) and 2(FF-PPC and SF-PPC) 



function getDataFigure1and2(sessionID,folderSourceString)

%sessionID: 1) arturo 2) wiggin 3) all ;

if ~exist('folderSourceString','var');    folderSourceString = 'E:/Mayo';     end
neuronType = 'All';
populationType = 'Stimulated';
oriChangeList = [2 3];
tapers=[2 3];

sessionNameString = getAttentionExperimentDetails;
if strcmp(sessionID,'all')
    fileNameStringList = cat(2,sessionNameString{1},sessionNameString{2});
elseif strcmp(sessionID,'arturo')
    fileNameStringList = sessionNameString{1};
elseif stcmp(sessionID,'wiggin')
    fileNameStringList = sessionNameString{2};
end

folderSave = fullfile(folderSourceString,'Data','dataForFig1and2');
makeDirectory(folderSave);
for i=1:length(fileNameStringList)
    
    disp(['Working on ' num2str(i) ' of ' num2str(length(fileNameStringList)) ' : ' fileNameStringList{i}])
    saveData(fileNameStringList{i},folderSourceString,folderSave,neuronType,populationType,oriChangeList,tapers);
    
end
end

function saveData(fileNameString,folderSourceString,folderSave,neuronType,populationType,oriChangeList,tapers)


attCueList = [{'H0V'} {'H1V'} {'HN'}];
folderName = fullfile(folderSourceString,'Data','segmentedData',fileNameString);
tpStr{1} = 'Baseline';      tpStr{2} = 'StimOnset';     tpStr{3} = 'TargetOnset';     
timePeriod{1} = [-0.25 0];      timePeriod{2} = [0.25 0.5];     timePeriod{3} = [-0.5 0];
binWidthMS = 10;
numCondition = length(attCueList);
electrodeArray = getGoodElectrodes(folderSourceString,fileNameString,neuronType,populationType);
[electrodePairsWithinHemisphere, electrodePairsAcrossHemispheres] = getGoodElectrodePairs(electrodeArray);
[~,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString);
goodList = [goodIndexList(1:2) goodIndexList(9)]; % Hit conditions for 0V, 1V and N

for i=1:length(tpStr)
    clear psdData psthData erpData FRData ffPPC sfPPC freqVals freqValsPPC timeVals xsFR
    fileNameSave = fullfile(folderSave,[fileNameString neuronType populationType 'Ori' num2str(oriChangeList(1)) num2str(oriChangeList(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '_' tpStr{i} '.mat']);
    fileNameSaveCoherence = fullfile(folderSave,[fileNameString neuronType populationType 'Ori' num2str(oriChangeList(1)) num2str(oriChangeList(2)) '_tapers' num2str(tapers(1)) num2str(tapers(2)) '_' tpStr{i} '_Coherence.mat']);
    if exist(fileNameSave,'file') && exist(fileNameSaveCoherence,'file')
        continue
    else
        for j=1:numCondition
            if i==1 || i==2
                lfpData = load(fullfile(folderName,[fileNameString attCueList{j} '_StimOnset_LFP']));
                spikeData = load(fullfile(folderName,[fileNameString attCueList{j} '_StimOnset_Spikes']));
            elseif i==3
                lfpData = load(fullfile(folderName,[fileNameString attCueList{j} '_TargetOnset_LFP']));
                spikeData = load(fullfile(folderName,[fileNameString attCueList{j} '_TargetOnset_Spikes']));
            end
            timeVals = lfpData.timeVals;
            timePos = intersect(find(timeVals>=timePeriod{i}(1)),find(timeVals<timePeriod{i}(2)));
            % Get Good Position
            x = orientationChangeDeg(goodList{j});
            goodPos = [];
            for m=1:length(oriChangeList)
                goodPos = cat(2,goodPos,find(x==uniqueOrientationChangeDeg(oriChangeList(m))));
            end
            goodPos = sort(goodPos);
            % MT analysis parameters for PSD
            Fs              = round(1/(timeVals(2)-timeVals(1)));
            params.tapers   = tapers;
            params.pad      = -1;
            params.Fs       = Fs;
            params.fpass    = [0 100];
            params.trialave = 1;
            
            for side=1:2
                eList = electrodeArray{side};
                ePairList = electrodePairsWithinHemisphere{side};
                for k=1:length(eList)
                    % PSD
                    [psdData{side}(j,k,:),freqVals] = mtspectrumc(squeeze(lfpData.segmentedLFPData(eList(k),goodPos,timePos))',params); %#ok<ASGLU>
                    % PSTH
                    [psthData{side}(j,k,:),xsFR] = getPSTH(spikeData.segmentedSpikeData(eList(k),goodPos),binWidthMS,[timeVals(1) timeVals(end)]);        %#ok<ASGLU>
                    FRData{side}(j,k) = mean(getSpikeCounts(spikeData.segmentedSpikeData(eList(k),goodPos),timePeriod{i}))/diff(timePeriod{i});
                end
                %ERP
                erpData{side}(j,:,:) = squeeze(mean(lfpData.segmentedLFPData(eList,goodPos,:),2));
                %FFPPC and SFPPC within hemisphere
                [ffPPC{side}(j,:,:),sfPPC{side}(j,:,:),freqValsPPC] = getPairWisePhaseConsistencyMeasures(lfpData.segmentedLFPData(:,goodPos,timePos),spikeData.segmentedSpikeData(:,goodPos),ePairList,tapers,timeVals,timePeriod{i}); %#ok<ASGLU>
            end
            % FFPPC and SFPPC across hemispheres
            [ffPPC{4}(j,:,:),sfPPC{4}(j,:,:)] = getPairWisePhaseConsistencyMeasures(lfpData.segmentedLFPData(:,goodPos,timePos),spikeData.segmentedSpikeData(:,goodPos),electrodePairsAcrossHemispheres,tapers,timeVals,timePeriod{i});
        end
        psdData{3} = combineData(psdData);
        psthData{3} = combineData(psthData);
        FRData{3} = combineData(FRData);
        erpData{3} = combineData(erpData);
        ffPPC{3} = combineData(ffPPC);
        sfPPC{3} = combineData(sfPPC);
              
        % Save file
        save(fileNameSave,'psdData','freqVals','psthData','FRData','xsFR','erpData','timeVals')
        save(fileNameSaveCoherence,'ffPPC','sfPPC','freqValsPPC')
    end
end
end

function electrodeArray = getGoodElectrodes(folderSourceString,fileNameString,neuronType,populationType)

suaCutoff = 3;
% Get sorting rating
sortRatings=sortRating(fileNameString,folderSourceString);
sua=intersect(find(sortRatings>0),find(sortRatings<=suaCutoff));
mua=find(sortRatings>suaCutoff);
all=union(sua,mua,'sorted');

[~,~,electrodeArrayPos]=electrodePositionOnGridMayo(1,fileNameString);

spkData{1}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H1V_StimOnset_Spikes']));
spkData{2}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H0V_StimOnset_Spikes']));

if strcmp(neuronType,'SUA')
    electrodeArray{1}=intersect(sua,electrodeArrayPos(:,8:13)); % Right Array
    electrodeArray{2}=intersect(sua,electrodeArrayPos(:,1:6)); % Left Array
    
elseif strcmp(neuronType,'MUA')
    electrodeArray{1}=intersect(mua,electrodeArrayPos(:,8:13));
    electrodeArray{2}=intersect(mua,electrodeArrayPos(:,1:6));
    
elseif strcmp(neuronType,'All')
    electrodeArray{1}=intersect(all,electrodeArrayPos(:,8:13));
    electrodeArray{2}=intersect(all,electrodeArrayPos(:,1:6));
end

if strcmp(populationType,'Stimulated')
    for j=1:2 %for each array
        elecList=electrodeArray{j}(:);
        for k=1:length(elecList)
            FRBl(j,k) =mean(getSpikeCounts(spkData{j}.segmentedSpikeData(elecList(k),:),[-0.25 0]))/diff([-0.25 0]); 
            FRSt(j,k) =mean(getSpikeCounts(spkData{j}.segmentedSpikeData(elecList(k),:),[0.25 0.5]))/diff([0.25 0.5]); 
        end
        delFR{j}=FRSt(j,:)-FRBl(j,:); 
        electrodeArray{j}=elecList(find(delFR{j}>5)); %#ok<FNDSB> % Threshold for choosing stimulated units is 5 spikes/s
    end
end
end

function [electrodePairsWithinHemisphere, electrodePairsAcrossHemispheres] = getGoodElectrodePairs(electrodeArray)

% get good Electrode pairs - within hemisphere
for j=1:2 % each electrode array
    electrodePairsWithinHemisphere{j} = combnk(electrodeArray{j},2); % Getting different electrode pairs for FFC and SFC using MATLAB combnk function
end

% get good Electrode pairs - across hemispheres
electrodePairsAcrossHemispheres = setdiff(combnk([electrodeArray{1}; electrodeArray{2}],2),[electrodePairsWithinHemisphere{1};electrodePairsWithinHemisphere{2}],'rows');

end

function [ffPPC,sfPPC,freqValsPPC] = getPairWisePhaseConsistencyMeasures(lfpData,spikeData,electrodePair,tapers,timeVals,timeRange)

% PPC- adapted from fieldtrip/connectivity/ft_connectivity_ppc.m computes 
% pairwise phase consistency  from a data-matrix containing a cross-spectral 
% density. This implements the method described in Vinck M, van Wingerden M, 
% Womelsdorf T, Fries P, Pennartz CM.
% The pairwise phase consistency: a bias-free measure of rhythmic neuronal
% synchronization. Vinck et al. Neuroimage. 2010

% Set up MT
Fs              = round(1/(timeVals(2)-timeVals(1)));
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 0;

% Field-Field Pairwise Phase Consistency
for i=1:size(electrodePair,1)
    clear lfp1 lfp2 input 
%    disp(['ffPPC ElectrodePair:',num2str(i)]);
    lfp1 = squeeze(lfpData(electrodePair(i,1),:,:));
    lfp2 = squeeze(lfpData(electrodePair(i,2),:,:));
    
    [~,~,S12_ffc,~,~,freqValsFFPPC]=coherencyc(lfp1',lfp2',params); %#ok<*AGROW>
    %ffPPC
    input = (S12_ffc./abs(S12_ffc))'; % normalize the cross-spectrum
    siz = size(input);
    n = siz(1);
    if n>1
        outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
        ffPPC(i,:)  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
    else
        error('computation of PPC requires >1 trial, please feed all trial dataset into computeCoherencyFromSpectrum program')
    end
end

% Spike-Field Pairwise Phase Consistency
nPairs = size(electrodePair,1);
if nPairs==nchoosek(length(unique(electrodePair)),2)     % for electrodePairsWithinHemisphere
    if isrow(unique(electrodePair))
        electrodePair = [electrodePair ; electrodePair(:,[2,1]) ; repmat(unique(electrodePair)',1,2)];   % it includes same electrode as pairs (d=0);
    else
        electrodePair = [electrodePair ; electrodePair(:,[2,1]) ; repmat(unique(electrodePair),1,2)];
    end
    if length(electrodePair)~=(2*nPairs)+length(unique(electrodePair))
        error('Number of within hemisphere electrode pairs are wrong')
    end
else
    electrodePair = [electrodePair ; electrodePair(:,[2,1])];   % for electrodePairsAcrossHemisphere
    if length(electrodePair)~= 2*nPairs
        error('Number of across hemisphere electrode pairs are wrong')
    end
end
for i=1:size(electrodePair,1)
    clear lfp spk input
    %    disp(['sfPPC ElectrodePair:',num2str(i)]);
    lfp = squeeze(lfpData(electrodePair(i,1),:,:));
    spk = convertSpikeTimes2Bins(spikeData(electrodePair(i,2),:,:),timeRange,1000/Fs);
    [~,~,S12_sfc,~,~,freqValsSFPPC]=coherencycpb(lfp',spk,params); %#ok<*AGROW>
    %sfPPC
    input = (S12_sfc./abs(S12_sfc))'; % normalize the cross-spectrum
    siz = size(input);
    n = siz(1);
    if n>1
        outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
        sfPPC(i,:)  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
    else
        error('computation of PPC requires >1 trial, please feed all trial dataset into computeCoherencyFromSpectrum program')
    end
end

% Sanity Check
if isequal(freqValsFFPPC,freqValsSFPPC)
    freqValsPPC = freqValsFFPPC;
else
    error('freqVals from FF-PPC & SF-PPC do not match!')
end

end

function sortRating=sortRating(fileNameString,folderSourceString)

[~,~,ratings]=xlsread(fullfile(folderSourceString,'Data','extractedData','SortRatings_CUonly'),fileNameString,'A1:A96');
if length(ratings)~=96
    error('electrode numbers do not match')
end

for i=1:96
    if ischar(ratings{i})
        ratings{i}=NaN;       % to eliminate electrodes rated as x, xx and blank.
    end
end
sortRating=cell2mat(ratings);
end
function combinedData = combineData(data)
if numel(size(data{1}))==3
    Data{1} = cat(2,data{1}(1,:,:),data{2}(2,:,:)); % Attend In - Valid [(R)H0V & (L)H1V]
    Data{2} = cat(2,data{1}(2,:,:),data{2}(1,:,:)); % Attend Out - Valid [(R)(H1V & (L)H0V)]
    Data{3} = cat(2,data{1}(3,:,:),data{2}(3,:,:)); % Neutral    [(R)N & (L)N]
    
elseif numel(size(data{1}))==2
    Data{1} = cat(2,data{1}(1,:),data{2}(2,:)); % Attend In - Valid [(R)H0V & (L)H1V]
    Data{2} = cat(2,data{1}(2,:),data{2}(1,:)); % Attend Out - Valid [(R)(H1V & (L)H0V)]
    Data{3} = cat(2,data{1}(3,:),data{2}(3,:)); % Neutral    [(R)N & (L)N]
end
combinedData = cat(1,Data{1},Data{2},Data{3});
end