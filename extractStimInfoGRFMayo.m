% This is copied from extractDigitalDataGRFLL. This reads all the
% stimulus information from the LL file

function extractStimInfoGRFMayo(fileNameString,folderSourceString,activeSide)

if ~exist('folderSourceString','var');   folderSourceString='F:\Projects\MayoProject\';  end

folderNameOut = fullfile(folderSourceString,'Data','segmentedData',fileNameString);

readDigitalCodesGRF(folderNameOut,activeSide); % writes stimResults
getGoodStimNumsGRF(folderNameOut,activeSide); % Good stimuli

end
% GRF Specific protocols
function stimResults = readDigitalCodesGRF(folderNameOut,activeSide)

% Load Lablib data structure
load(fullfile(folderNameOut,'LL.mat'));

%%%%%%%%%%%%%%%%% Get info from LL to construct stimResults %%%%%%%%%%%%%%%
if activeSide==0 % Map0
    validMap = find(LL.stimType1>0);
    aziLL = LL.azimuthDeg1(validMap);
    eleLL = LL.elevationDeg1(validMap);
    sigmaLL = LL.sigmaDeg1(validMap);
    
    if isfield(LL,'radiusDeg1')
        radiusExists = 1;
        radiusLL = LL.radiusDeg1(validMap);
    else
        radiusExists = 0;
    end
    sfLL = LL.spatialFreqCPD1(validMap);
    oriLL = LL.orientationDeg1(validMap);
    conLL = LL.contrastPC1(validMap);  
    timeLL = LL.time1(validMap)/1000;
    taskType = LL.stimType1(validMap);
    
elseif activeSide==1 % Map2
    
    validMap = find(LL.stimType2>0);
    aziLL = LL.azimuthDeg2(validMap);
    eleLL = LL.elevationDeg2(validMap);
    sigmaLL = LL.sigmaDeg2(validMap);
    if isfield(LL,'radiusDeg2')
        radiusExists = 1;
        radiusLL = LL.radiusDeg2(validMap);
    else
        radiusExists = 0;
    end
    sfLL = LL.spatialFreqCPD2(validMap);
    oriLL = LL.orientationDeg2(validMap);
    conLL = LL.contrastPC2(validMap);
    timeLL = LL.time2(validMap)/1000;
    taskType = LL.stimType2(validMap);
end

stimResults.azimuth = aziLL;
stimResults.elevation = eleLL;
stimResults.contrast = conLL;
if radiusExists
    stimResults.radius = radiusLL;
end
stimResults.sigma = sigmaLL;
stimResults.orientation = oriLL;
stimResults.spatialFrequency = sfLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTrials = length(LL.eotCode);
trialStartTimesLL = LL.startTime;
instructionTrials = LL.instructTrial;
catchTrials = LL.catchTrial;
trialCertify = LL.trialCertify;
eotCodes = LL.eotCode;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numStims  = getStimPosPerTrial(trialStartTimesLL,timeLL);
stimResults.side = activeSide;

posTask = 0;
pos=0;
for i=1:numTrials 
    
    taskTypeThisTrial = taskType(posTask+1:posTask+numStims(i));
    if (numStims(i)>0)
        stimTimesFromTrialStart = timeLL(pos+1:pos+numStims(i)) - trialStartTimesLL(i);
        stimResults.time(pos+1:pos+numStims(i)) = stimTimesFromTrialStart; % times relative to trialStart
        stimResults.type(pos+1:pos+numStims(i)) = taskTypeThisTrial;
        stimResults.trialNumber(pos+1:pos+numStims(i)) = i;
        stimResults.stimPosition(pos+1:pos+numStims(i)) = 1:numStims(i);
        
        stimResults.instructionTrials(pos+1:pos+numStims(i)) = instructionTrials(i); %always zero
        stimResults.catch(pos+1:pos+numStims(i)) = catchTrials(i);
        stimResults.eotCodes(pos+1:pos+numStims(i)) = eotCodes(i);
        stimResults.trialCertify(pos+1:pos+numStims(i)) = trialCertify(i);
        pos = pos+numStims(i);
    end
    posTask = posTask+numStims(i);
end

% Save in folderNameOut
save(fullfile(folderNameOut,['stimResults' num2str(activeSide) '.mat']),'stimResults');

end
function [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderNameOut,activeSide)

load(fullfile(folderNameOut,['stimResults' num2str(activeSide) '.mat']));

totalStims = length(stimResults.eotCodes);
disp(['Number of trials: ' num2str(max(stimResults.trialNumber))]);
disp(['Number of stimuli: ' num2str(totalStims)]);

% exclude uncertified trials, catch trials and instruction trials
tc = find(stimResults.trialCertify==1);
it = find(stimResults.instructionTrials==1);
ct = find(stimResults.catch==1);

badStimNums = [it tc ct];

%eottypes
% 0 - correct, 1 - wrong, 2-failed, 3-broke, 4-ignored, 5-False
% Alarm/quit, 6 - distracted, 7 - force quit
%disp('Analysing correct, wrong and failed trials');
%badEOTs = find(stimResults.eotCodes>2); 
%disp('Analysing correct and wrong trials')
%badEOTs = find(stimResults.eotCodes>1); 
disp('Analysing only correct trials')
badEOTs = find(stimResults.eotCodes>0); 
badStimNums = [badStimNums badEOTs];

goodStimNums = setdiff(1:totalStims,unique(badStimNums));

% stim types
% 0 - Null, 1 - valid, 2 - target, 3 - frontpadding, 4 - backpadding

disp('Only taking valid stims ');
validStims = find(stimResults.type==1);
goodStimNums = intersect(goodStimNums,validStims);

% %%%%%%%%%%%%%% Remove bad stimuli after target %%%%%%%%%%%%%%%%%%%%
% 
% clear trialNums stimPos
% trialNums = stimResults.trialNumber(goodStimNums);
% stimPos   = stimResults.stimPosition(goodStimNums);
% 
% % Get the target positions of the trialNums
% clear goodTrials
% goodTrials = unique(trialNums);
% 
% clear targetPos
% for i=1:length(goodTrials)
%     allStimWithThisTrialNum = find(stimResults.trialNumber==goodTrials(i));
%     
%     if sum(stimResults.catch(allStimWithThisTrialNum))>0        % catch trials
%         targetPos(trialNums==goodTrials(i)) = inf; %#ok<*AGROW>
%     else
%         targetPos(trialNums==goodTrials(i)) = find(stimResults.type(allStimWithThisTrialNum)==2);
%     end
% end
% 
% validStimuliAfterTarget = find(stimPos>targetPos);
% if ~isempty(validStimuliAfterTarget)
%     disp([num2str(length(validStimuliAfterTarget)) ' out of ' num2str(length(goodStimNums)) ' stimuli after target']);
%     save(fullfile(folderNameOut,'validStimAfterTarget.mat'),'validStimuliAfterTarget');
% end
% 
% goodStimNums(validStimuliAfterTarget)=[];
disp(['Number of good stimuli: ' num2str(length(goodStimNums))]);
goodStimTimes = stimResults.time(goodStimNums);
save(fullfile(folderNameOut,['goodStimNums' num2str(activeSide) '.mat']),'goodStimNums','goodStimTimes');

end
function [numStim,stimOnPos] = getStimPosPerTrial(trialStartTimes, stimStartTimes)

numTrials = length(trialStartTimes);

stimOnPos = cell(1,numTrials);
numStim   = zeros(1,numTrials);

for i=1:numTrials-1
    stimOnPos{i} = intersect(find(stimStartTimes>=trialStartTimes(i)),find(stimStartTimes<trialStartTimes(i+1)));
    numStim(i) = length(stimOnPos{i});
end
stimOnPos{numTrials} = find(stimStartTimes>=trialStartTimes(numTrials));
numStim(numTrials) = length(stimOnPos{numTrials});
end