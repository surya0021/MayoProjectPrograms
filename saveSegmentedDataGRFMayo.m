% This program takes the fileNameString and saves the segmented data in
% appropriate folders

function saveSegmentedDataGRFMayo(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='F:\Projects\MayoProject\';       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 2000;
timePeriodS = [0 0.4]; % Save this interval around each stimulus onset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderNameIn = fullfile(folderSourceString,'Data','extractedData');
folderNameOut = fullfile(folderSourceString,'Data','segmentedData',fileNameString);

fileNameLFP = strcat(fileNameString, '_extractedTrialsLFP_RFmap');
lfpData = load(fullfile(folderNameIn,fileNameLFP));
lfpData = lfpData.(fileNameLFP); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
display([fileNameLFP ' file loaded.....']);

s0=load(fullfile(folderNameOut,'stimResults0.mat'));
numStimuli = length(s0.stimResults.time);

%%%%%%%%%%%%%%%%%%%%%%%% Find appropriate indices %%%%%%%%%%%%%%%%%%%%%%%%%

numElectrodes = size(lfpData,2);
numTimePos=round(Fs*diff(timePeriodS));

segmentedLFPData = zeros(numElectrodes,numStimuli,numTimePos);

for i=1:numStimuli
    tmpLFP = cell2mat(lfpData(s0.stimResults.trialNumber(i),:)');
    startPos = round(Fs*s0.stimResults.time(i));
    
    segmentedLFPData(:,i,:) = tmpLFP(:,startPos+ (1:numTimePos));
end
timeVals = timePeriodS(1):1/Fs:timePeriodS(2)-1/Fs; %#ok<NASGU>
fileNameSave = fullfile(folderNameOut,'lfpData');
save(fileNameSave,'timeVals','segmentedLFPData');
end