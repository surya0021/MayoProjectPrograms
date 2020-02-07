%This function saves the neutral trials which are divided further depending
%on target location
%Here the format is {trialOutcome}{targetLocation}{cueType} used instead of
%previous format {trialOutcome}{cueLocation}{cueType} c.f.
%savesegmentedDataMayo.m
%{trialOutcome}= Hit(H) or Miss(M);  {targetLocation}= Left(0) or Right(1);
%{cueType}= Neutral(N)

function saveSegmentedNeutralDataMayo(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='E:\Mayo';       end
neutralTrialFlag=1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
saveStringConditionList = [{'H0N'} {'H1N'} {'M0N'} {'M1N'}]; % Data must be stored in this order
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
if isfield(lfpData,fileNameLFP)
    lfpData = lfpData.(fileNameLFP); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
else
    lfpData = lfpData.extractedTrialsLFP;
end
display([fileNameLFP ' file loaded.....']);

fileNameSpikes = strcat(fileNameString, '_Spikes');
spikeData = load(fullfile(folderNameIn,fileNameSpikes));
spikeData = spikeData.(fileNameSpikes); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
display([fileNameSpikes ' file loaded.....']);

%%%%%%%%%%%%%%%%%%%%%%%% Find appropriate indices %%%%%%%%%%%%%%%%%%%%%%%%%
goodIndexList = getGoodIndices(CDS,DAT,[],neutralTrialFlag); % Gets Indices for the 12 categories
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

for i=1:length(saveStringConditionList)
    disp(['Working on condition ' saveStringConditionList{i}]);
    
    for j=1:numTimeSegments % For each of the time periods
        clear segmentedLFPData segmentedSpikeData
                
        fileNameSaveLFP = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList{j} '_LFP']);
        fileNameSaveSpikes = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList{j} '_Spikes']);
        
        numTrials = length(goodIndexList{8+i});  % GoodIndexList for Neutral Trials starts from the 9th cell 9:H0N , 10:H1N , 11:M0N, 12:M1N
        segmentedLFPData = zeros(numElectrodes,numTrials,numTimePos(j));
        segmentedSpikeData = cell(numElectrodes,numTrials);
        
        for k=1:numTrials
            t=goodIndexList{8+i}(k); % Trial Number of interest
 
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