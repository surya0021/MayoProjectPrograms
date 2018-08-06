% This program is modified from saveLLData. Saves data from the DAT file
% provided by Patrick Mayo.

% All data is saved in segmentedData/fileNameString.

function LL = saveLLDataGRFMayo(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='F:\Projects\MayoProject\';  end

folderNameIn = fullfile(folderSourceString,'Data','extractedData');
fileNameDAT = strcat(fileNameString, '_DAT');
DAT = load(fullfile(folderNameIn,fileNameDAT));
DAT = DAT.(fileNameDAT);
DAT2 = convertDATFile(DAT);

folderNameOut = fullfile(folderSourceString,'Data','segmentedData',fileNameString);
makeDirectory(folderNameOut);

% Stimulus properties
numTrials = length(DAT2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allStimulusIndex = [];
allStimulusOnTimes = [];
gaborIndex = [];
stimType = [];
azimuthDeg = [];
elevationDeg = [];
sigmaDeg = [];
radiusDeg = [];
spatialFreqCPD =[];
orientationDeg = [];
contrastPC = [];

eotCode=[];
myEotCode=[];
endTime=[];
startTime=[];
instructTrial=[];
catchTrial=[];
trialCertify=[];

for i=1:numTrials
    disp(['trial ' num2str(i) ' of ' num2str(numTrials)]);
    clear trials
    trials = DAT2{i};
    
    if isfield(trials,'trial')
        instructTrial = [instructTrial [trials.trial.data.instructTrial]]; %#ok<*AGROW>
        catchTrial    = [catchTrial [trials.trial.data.catchTrial]];
    end
    
    if isfield(trials,'trialCertify')
        trialCertify = [trialCertify [trials.trialCertify.data]];
    end
    
    if isfield(trials,'trialEnd')
        eotCode = [eotCode [trials.trialEnd.data]];
        endTime = [endTime [trials.trialEnd.timeMS]];
    end
    
    if isfield(trials,'myTrialEnd')
        myEotCode = [myEotCode [trials.myTrialEnd.data]];
    end
    
    if isfield(trials,'trialStart')
        startTime = [startTime [trials.trialStart.timeMS]];
    end
    
    if isfield(trials,'stimulusOn')
        allStimulusIndex = [allStimulusIndex [trials.stimulusOn.data]'];
    end
    
    if isfield(trials,'stimulusOnTime')
        allStimulusOnTimes = [allStimulusOnTimes [trials.stimulusOnTime.timeMS]'];
    end
    
    if isfield(trials,'stimDesc')
        gaborIndex = [gaborIndex [trials.stimDesc.data.gaborIndex]];
        stimType = [stimType [trials.stimDesc.data.stimType]];
        azimuthDeg = [azimuthDeg [trials.stimDesc.data.azimuthDeg]];
        elevationDeg = [elevationDeg [trials.stimDesc.data.elevationDeg]];
        sigmaDeg = [sigmaDeg [trials.stimDesc.data.sigmaDeg]];
        if isfield(trials.stimDesc.data,'radiusDeg')
            radiusExists=1;
            radiusDeg = [radiusDeg [trials.stimDesc.data.radiusDeg]];
        else
            radiusExists=0;
        end
        spatialFreqCPD = [spatialFreqCPD [trials.stimDesc.data.spatialFreqCPD]];
        orientationDeg = [orientationDeg [trials.stimDesc.data.directionDeg]];
        contrastPC = [contrastPC [trials.stimDesc.data.contrastPC]];
    end
end

% Sort stim properties by stimType
numGabors = length(unique(gaborIndex));
for i=1:numGabors
    gaborIndexFromStimulusOn{i} = find(allStimulusIndex==i-1);
    gaborIndexFromStimDesc{i} = find(gaborIndex==i-1);
end

if isequal(gaborIndexFromStimDesc,gaborIndexFromStimulusOn)
    for i=1:numGabors
        LL.(['time' num2str(i-1)]) = allStimulusOnTimes(gaborIndexFromStimulusOn{i});
        LL.(['stimType' num2str(i-1)]) = stimType(gaborIndexFromStimulusOn{i});
        LL.(['azimuthDeg' num2str(i-1)]) = azimuthDeg(gaborIndexFromStimulusOn{i});
        LL.(['elevationDeg' num2str(i-1)]) = elevationDeg(gaborIndexFromStimulusOn{i});
        LL.(['sigmaDeg' num2str(i-1)]) = sigmaDeg(gaborIndexFromStimulusOn{i});
        if radiusExists
            LL.(['radiusDeg' num2str(i-1)]) = radiusDeg(gaborIndexFromStimulusOn{i});
        end
        LL.(['spatialFreqCPD' num2str(i-1)]) = spatialFreqCPD(gaborIndexFromStimulusOn{i});
        LL.(['orientationDeg' num2str(i-1)]) = orientationDeg(gaborIndexFromStimulusOn{i});
        LL.(['contrastPC' num2str(i-1)]) = contrastPC(gaborIndexFromStimulusOn{i});
    end
else
    error('Gabor indices from stimuluOn and stimDesc do not match!!');
end

LL.eotCode = eotCode;
LL.myEotCode = myEotCode;
LL.startTime = startTime/1000; % in seconds
LL.endTime = endTime/1000; % in seconds
LL.instructTrial = instructTrial;
LL.catchTrial = catchTrial;
LL.trialCertify = trialCertify;

save(fullfile(folderNameOut,'LL.mat'),'LL');
end

function DAT2 = convertDATFile(DAT)
i=1; continueFlag=1;
while(continueFlag)
    tnum = (['t', num2str(i)]);
    if isfield(DAT,tnum)
        DAT2{i} = DAT.(tnum); % we keep only the data about the kept trials
        i=i+1;
    else
        continueFlag=0;
    end
end
end