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