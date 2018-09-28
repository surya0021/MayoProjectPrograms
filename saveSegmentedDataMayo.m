% This program takes the fileNameString and saves the segmented data in
% appropriate folders

% fileNameString: string indicating the fileName.
% We assume that the following files are available at
% folderSourceString\Data\extractedData
% 1. fileNameString_Codes
% 2. fileNameString_DAT
% 3. fileNameString_extractedTrialsLFP

% Followng Patrick's nomenclature, we save data in the format
% {trialOutcome}{CueType}{CueLocation}.
% {trialOutcome}: Hit (H) or Miss (M)
% {CueType}: Valid (V), Invalid (I) or Neutral (N)
% {CueLocation}: In case the CueType is Valid or Invalid, we also have the
% cue location (0: left or 1: right).
% So there are the following 10 file types: H0V,H1V,H0I,H1I,M0V,M1V,M0I,M1I,HN,MN

% Further, we save data around two events of interest: first stimulus and the target.

function saveSegmentedDataMayo(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
saveStringConditionList = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'HN'} {'MN'}]; % Data must be stored in this order
Fs = 2000;
timePeriodS{1} = [-0.25 0.5]; % Save this interval around the first stimulus onset 
timePeriodS{2} = [-0.5 0.1]; % Save this interval around the target onset
saveStringTimePeriodList = [{'_StimOnset'} {'_TargetOnset'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderNameIn = fullfile(folderSourceString,'Data','extractedData');
folderNameOut = fullfile(folderSourceString,'Data','segmentedData',fileNameString);
makeDirectory(folderNameOut);

fileNameCDS = [fileNameString '_Codes'];
CDS = load(fullfile(folderNameIn,fileNameCDS));
CDS = CDS.(fileNameCDS);
display([fileNameCDS ' loaded....']);

fileNameDAT = strcat(fileNameString, '_DAT');
DAT = load(fullfile(folderNameIn,fileNameDAT));
DAT = DAT.(fileNameDAT);
display([fileNameDAT ' loaded....']);

fileNameLFP = strcat(fileNameString, '_extractedTrialsLFP');
lfpData = load(fullfile(folderNameIn,fileNameLFP));
lfpData = lfpData.(fileNameLFP); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
display([fileNameLFP ' file loaded.....']);

fileNameSpikes = strcat(fileNameString, '_Spikes');
spikeData = load(fullfile(folderNameIn,fileNameSpikes));
spikeData = spikeData.(fileNameSpikes); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
display([fileNameSpikes ' file loaded.....']);

%%%%%%%%%%%%%%%%%%%%%%%% Find appropriate indices %%%%%%%%%%%%%%%%%%%%%%%%%
goodIndexList = getGoodIndices(CDS,DAT); % Get Indices for the 10 categories
trialStartTimeS = cell2mat(cellfun(@(x) x{1}, CDS(:,1),'UniformOutput',false ));
stimulusOnTimeS = cellfun(@(x) x{5}, CDS(:,1),'UniformOutput',false); % still under cell format. dimension: trials x 1 (in the first and only columnm there are cells containing all the trial start times for the corresponding trial
%saccadeTimeS = cell2mat(cellfun(@(x) x{7}, CDS(:,1),'UniformOutput',false));

[~,~,targetOnTimeMS,~,~] = getInfoDATFile(DAT);

numElectrodes = size(lfpData,2);
numTimeSegments = length(timePeriodS);
numTimePos=zeros(1,numTimeSegments);
for i=1:numTimeSegments
    numTimePos(i) = round(Fs*diff(timePeriodS{i}));
end

for i=1:length(goodIndexList) % For each of the 10 conditions
    disp(['Working on condition ' saveStringConditionList{i}]);
    
    for j=1:numTimeSegments % For each of the time periods
        clear segmentedLFPData segmentedSpikeData
    
        fileNameSaveLFP = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList{j} '_LFP']);
        fileNameSaveSpikes = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList{j} '_Spikes']);
        numTrials = length(goodIndexList{i});
        segmentedLFPData = zeros(numElectrodes,numTrials,numTimePos(j));
        segmentedSpikeData = cell(numElectrodes,numTrials);
        
        for k=1:numTrials
            t=goodIndexList{i}(k); % Trial Number of interest
 
            if j==1
                startPosLFP = round(Fs*(stimulusOnTimeS{t}(1)-trialStartTimeS(t) + timePeriodS{j}(1)));
                startPosSpikes = stimulusOnTimeS{t}(1);
            elseif j==2
                startPosLFP = round(Fs*(targetOnTimeMS(t)/1000 + stimulusOnTimeS{t}(1)-trialStartTimeS(t) + timePeriodS{j}(1)));
                startPosSpikes = targetOnTimeMS(t)/1000 + stimulusOnTimeS{t}(1);
            end
        
            for n=1:numElectrodes
                segmentedLFPData(n,k,:) = lfpData{t,n}(startPosLFP+(1:numTimePos(j)));
                spikeDataThisTrial = spikeData{n,1,t};
                if ~isempty(spikeDataThisTrial)
                    spikeTimesThisTrial = spikeDataThisTrial - startPosSpikes; % times relative to stimulus or target onset
                    usefulSpikePos = intersect(find(spikeTimesThisTrial>=timePeriodS{j}(1)),find(spikeTimesThisTrial<timePeriodS{j}(2)));
                    segmentedSpikeData{n,k} = spikeTimesThisTrial(usefulSpikePos);
                end
            end
        end
        timeVals = timePeriodS{j}(1):1/Fs:timePeriodS{j}(2)-1/Fs; %#ok<NASGU>
        save(fileNameSaveLFP,'timeVals','segmentedLFPData');
        save(fileNameSaveSpikes,'segmentedSpikeData');
    end
end
end

function goodIndexList = getGoodIndices(CDS,DAT)

[DAT2,CueOnCode] = getInfoDATFile(DAT);
DAT2 = DAT2(1:length(CDS)); % DAT file sometimes has more entries than the CDS file. The last block for attLoc1 has not been used for analysis.
CueOnCode = CueOnCode(1:length(CDS));

%trialEndCode = cell2mat(cellfun(@(x) x.trialEnd.data, DAT2,'UniformOutput',false));
trialEndCode = cell2mat(cellfun(@(x) x{8},CDS(:,2),'UniformOutput',false))';

isValidTrial = cell2mat(cellfun(@(x) x.trial.data.validTrial,DAT2,'UniformOutput',false));
% isInstructTrial = cell2mat(cellfun(@(x) x.trial.data.instructTrial,DAT2,'UniformOutput',false));
isInstructTrial = (CueOnCode>0); % Sometimes a cue is generated even if the trial is not an instruction trial. Therefore, any trial which contains the cueOn field is considered an instruction trial
attendLoc = cell2mat(cellfun(@(x) x.trial.data.attendLoc,DAT2,'UniformOutput',false));
isCatchTrial = cell2mat(cellfun(@(x) x.trial.data.catchTrial,DAT2,'UniformOutput',false));

% Hit Trials
hitIndices = (trialEndCode==0) & (isInstructTrial==0) & (isCatchTrial==0);
goodIndexList{1} = find(hitIndices & (isValidTrial==1) & (attendLoc==0)); % H0V
goodIndexList{2} = find(hitIndices & (isValidTrial==1) & (attendLoc==1)); % H1V
goodIndexList{3} = find(hitIndices & (isValidTrial==0) & (attendLoc==0)); % H0I
goodIndexList{4} = find(hitIndices & (isValidTrial==0) & (attendLoc==1)); % H1I

% Miss Trials
missIndices = (trialEndCode==2) & (isInstructTrial==0) & (isCatchTrial==0);
goodIndexList{5} = find(missIndices & (isValidTrial==1) & (attendLoc==0)); % M0V
goodIndexList{6} = find(missIndices & (isValidTrial==1) & (attendLoc==1)); % M1V
goodIndexList{7} = find(missIndices & (isValidTrial==0) & (attendLoc==0)); % M0I
goodIndexList{8} = find(missIndices & (isValidTrial==0) & (attendLoc==1)); % M1I

goodIndexList{9} = find(hitIndices & (isValidTrial==1) & (attendLoc==2)); % HN
goodIndexList{10} = find(missIndices & (isValidTrial==1) & (attendLoc==2)); % MN

% Doing the same thing using Patrick's code
[~, indHitLoc0, indHitLoc1,indHitNeutral] = getTrialTypes (CDS,1,0,0);
catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indHitLoc0 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indHitLoc0 ); % logical index of valid trials for Hit Loc 0 trials
goodIndexList2{1} = (indHitLoc0(validList & ~catchList))'; % HOV
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indHitLoc0 ); % logical index of invalid trials for Hit Loc 0 trials
goodIndexList2{3} = (indHitLoc0(invalidList & ~catchList))'; % HOI

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indHitLoc1 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indHitLoc1 ); % logical index of valid trials for Hit Loc 1 trials
goodIndexList2{2} = (indHitLoc1(validList & ~catchList))'; % H1V
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indHitLoc1 ); % logical index of invalid trials for Hit Loc 1 trials
goodIndexList2{4} = (indHitLoc1(invalidList & ~catchList))'; % H1I

[~, indMissLoc0, indMissLoc1,indMissNeutral] = getTrialTypes (CDS,0,1,0);
catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indMissLoc0 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indMissLoc0 ); % logical index of valid trials for Miss Loc 0 trials
goodIndexList2{5} = (indMissLoc0(validList & ~catchList))'; % MOV
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indMissLoc0 ); % logical index of invalid trials for Miss Loc 0 trials
goodIndexList2{7} = (indMissLoc0(invalidList & ~catchList))'; % MOI

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indMissLoc1 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indMissLoc1 ); % logical index of valid trials for Miss Loc 1 trials
goodIndexList2{6} = (indMissLoc1(validList & ~catchList))'; % M1V
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indMissLoc1 ); % logical index of invalid trials for Miss Loc 1 trials
goodIndexList2{8} = (indMissLoc1(invalidList & ~catchList))'; % M1I

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indHitNeutral ); % logical index of catch trials with no stimulus change
goodIndexList2{9} = (indHitNeutral(~catchList))'; % HN

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indMissNeutral ); % logical index of catch trials with no stimulus change
goodIndexList2{10} = (indMissNeutral(~catchList))'; % MN

if ~isequal(goodIndexList,goodIndexList2)
    error('Index Lists do not match...');
end
end