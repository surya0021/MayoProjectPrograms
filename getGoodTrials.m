% This function gives the trial numbers which have the required trial conditions and orientation
% change.
% changeIndex: [1:6]; 1-least orientation change, 6-highest orientation
% change of the respective sessions.
% goodPos: Trial numbers as of in segmented data.
% goodTrials: Trial numbers as of in DAT file.



function [goodPos,goodTrials]= getGoodTrials(fileNameString,conditionString,changeIndex)

if ~exist('changeIndex','var'); changeIndex=[2,3]; end % orientation change considered in Mayo and Maunsell(2016) for analyzing the data.

folderSourceString='E:\Mayo'; folderNameIn = fullfile(folderSourceString,'Data','extractedData');
conditionStringList = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'HN'} {'MN'}];
for i=1:length(conditionStringList)
    if strcmp(conditionStringList{i},conditionString)
        conditionIndex=i;
    end
end

fileNameCDS = [fileNameString '_Codes'];
        CDS = load(fullfile(folderNameIn,fileNameCDS));
        CDS = CDS.(fileNameCDS);
fileNameDAT = strcat(fileNameString, '_DAT');
        DAT = load(fullfile(folderNameIn,fileNameDAT));
        DAT = DAT.(fileNameDAT);
     

goodIndexList=getGoodIndex(CDS,DAT);

DAT2=getInfoDATFile(DAT);
oriChangeVals=[];
for i=1:length(DAT2);
oriChangeVals=[oriChangeVals [DAT2{i}.trial.data.change]]; %#ok<AGROW>
end
oriChange=unique(oriChangeVals,'sorted');

oriChangeTrials=find(ismember(oriChangeVals,oriChange(changeIndex)));
goodTrials=intersect(goodIndexList{conditionIndex},oriChangeTrials);
goodPos=find(ismember(goodIndexList{conditionIndex},goodTrials));
     

end


function goodIndexList = getGoodIndex(CDS,DAT)

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


