function saveDataForLDAMayo(fileNameString,folderSourceString,neutralTrialFlag)

if ~exist('folderSourceString','var') || isempty(folderSourceString); folderSourceString='E:/Mayo';  end;
if ~exist('neutralTrialFlag','var');    neutralTrialFlag=0;     end

tapers=[2 3];   %tapers=[TW K]
tpStr='TargetOnset_500ms';
populationType='Stimulated';    % populationType (1)'All': all good electrodes of a session     (2)'Stimulated' : gives the electrodes which are responsive (DelFR > 5Hz).
timeRange=[-0.5 0];

folderSave = fullfile(folderSourceString,'Data','savedDataForLDA');
makeDirectory(folderSave);
if neutralTrialFlag
    fileToSave = fullfile(folderSave,[fileNameString populationType tpStr '_Tapers_' num2str(tapers(2)) '_SplitNeutral.mat']);
else
    fileToSave = fullfile(folderSave,[fileNameString populationType tpStr '_Tapers_' num2str(tapers(2)) '.mat']);
end
% Load Data

folderName=fullfile(folderSourceString,'Data','segmentedData',fileNameString);
if neutralTrialFlag
    attCueList={'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'};
else
    attCueList={'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','HN','MN'};
end
for i=1:length(attCueList)
    lfpData{i}=load(fullfile(folderName,[fileNameString attCueList{i} '_TargetOnset_LFP']));
    spikeData{i}=load(fullfile(folderName,[fileNameString attCueList{i} '_TargetOnset_Spikes'])); 
end

timeVals=lfpData{1}.timeVals;
Fs = round(1/(timeVals(2)-timeVals(1)));
pos = find(timeVals>=timeRange(1),1)+ (0:diff(timeRange)*Fs-1);

[eList,RightArrayPos,LeftArrayPos]=getGoodElectrodes(fileNameString,folderSourceString,populationType); %#ok<ASGLU>
[~,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString,neutralTrialFlag);

freqVals = 0:1/(diff(timeRange)):Fs-1/(diff(timeRange));
alphaRangeHz = [8 12]; gammaRangeHz = [40 80]; lineNoiseFreqHz= 60;
alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
lineNoisePos = find(freqVals==lineNoiseFreqHz);
gammaPos = setdiff(intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<=gammaRangeHz(2))),lineNoisePos);

% Parameters for MT

params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 0;

for ori=1:length(uniqueOrientationChangeDeg)
    disp(['Orientation ' num2str(ori) ':'])
    for cond=1:length(attCueList)
        clear oriChangeDegThisCondition selectedPos -regexp ^analyzed
        oriChangeDegThisCondition=orientationChangeDeg(goodIndexList{cond});
        selectedPos=find(oriChangeDegThisCondition==uniqueOrientationChangeDeg(ori));
        if isempty(selectedPos)
            firingRate{ori}{cond}=[];
            PSDData{ori}{cond}=[];
            alphaData{ori}{cond}=[];
            gammaData{ori}{cond}=[];
            disp([attCueList{cond} ' Stim Repeats ' num2str(0)]);
            continue
        end
        analyzedLFPData=lfpData{cond}.segmentedLFPData(eList,selectedPos,pos);
        analyzedSpikeData=spikeData{cond}.segmentedSpikeData(eList,selectedPos);
        disp([attCueList{cond} ' Stim Repeats ' num2str(size(analyzedLFPData,2))]);
        for k=1:length(eList)
            firingRate{ori}{cond}(k,:)=getSpikeCounts(analyzedSpikeData(k,:),timeRange)/diff(timeRange);
            [PSDData{ori}{cond}(k,:,:),freqValsMT]= mtspectrumc(squeeze(analyzedLFPData(k,:,:))',params);
        end
        if ~isequal(unique(diff(freqVals)),unique(diff(freqValsMT)))
            error('Frequency values does not match')
        end
        alphaData{ori}{cond}=squeeze(sum(PSDData{ori}{cond}(:,alphaPos,:),2));
        gammaData{ori}{cond}=squeeze(sum(PSDData{ori}{cond}(:,gammaPos,:),2));
    end
end
save(fileToSave,'firingRate','freqValsMT','PSDData','alphaData','gammaData','RightArrayPos','LeftArrayPos','attCueList','uniqueOrientationChangeDeg');
end



function [eList,RightArrayPos,LeftArrayPos]=getGoodElectrodes(fileNameString,folderSourceString,populationType)

suaCutoff = 3;
% Get sorting rating
sortRatings=sortRating(fileNameString,folderSourceString);
sua=intersect(find(sortRatings>0),find(sortRatings<=suaCutoff));
mua=find(sortRatings>suaCutoff);
all=union(sua,mua,'sorted');
[~,~,electrodeArrayPos]=electrodePositionOnGridMayo(1,fileNameString);
electrodeArray{1}=intersect(all,electrodeArrayPos(:,8:13)); % Right array
electrodeArray{2}=intersect(all,electrodeArrayPos(:,1:6)); % Left array

if strcmp(populationType,'Stimulated')
    spkData{1}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H1V_StimOnset_Spikes']));
    spkData{2}=load(fullfile(folderSourceString,'Data','segmentedData',fileNameString,[fileNameString 'H0V_StimOnset_Spikes']));
    for j=1:2 %for each array
        elecList=electrodeArray{j}(:);
        for k=1:length(elecList)
            FRBl(j,k) =mean(getSpikeCounts(spkData{j}.segmentedSpikeData(elecList(k),:),[-0.25 0]))/diff([-0.25 0]);
            FRSt(j,k) =mean(getSpikeCounts(spkData{j}.segmentedSpikeData(elecList(k),:),[0.25 0.5]))/diff([0.25 0.5]);
        end
        delFR{j}=FRSt(j,:)-FRBl(j,:);
        electrodeArray{j}=elecList(find(delFR{j}>5)); %#ok<FNDSB> % Threshold for choosing stimulated units is 5 spikes/s
    end
    eList=sort(cat(1,electrodeArray{:}));
elseif strcmp(populationType,'All');
    eList=all;
end
[~,RightArrayPos]=intersect(eList,electrodeArrayPos(:,8:13));
[~,LeftArrayPos]=intersect(eList,electrodeArrayPos(:,1:6));
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