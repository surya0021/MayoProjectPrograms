function saveDataForLDAMayo(fileNameString,folderSourceString)



if ~exist('folderSourceString','var');  folderSourceString='E:/Mayo';  end;
tapers=[2 3];   %tapers=[TW K]
tpStr='TargetOnset_500ms';
populationType='All';
timeRange=[-0.5 0];

folderSave = fullfile(folderSourceString,'Data','savedDataForLDA');
makeDirectory(folderSave);
fileToSave = fullfile(folderSave,[fileNameString populationType tpStr '_Tapers_' num2str(tapers(2)) '.mat']);

% Load Data

folderName=fullfile(folderSourceString,'Data','segmentedData',fileNameString);
attCueList={'H0V','H1V','H0I','H1I','HN','M0V','M1V','M0I','M1I','MN'};
for i=1:10
    lfpData{i}=load(fullfile(folderName,[fileNameString attCueList{i} '_TargetOnset_LFP']));
    spikeData{i}=load(fullfile(folderName,[fileNameString attCueList{i} '_TargetOnset_Spikes'])); 
end

timeVals=lfpData{1}.timeVals;
Fs = round(1/(timeVals(2)-timeVals(1)));
pos = find(timeVals>=timeRange(1),1)+ (0:diff(timeRange)*Fs-1);

[eList,RightArrayPos,LeftArrayPos]=getGoodElectrodes(fileNameString,folderSourceString); %#ok<ASGLU>
[~,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString);
goodList = [goodIndexList(1:4) goodIndexList(9) goodIndexList(5:8) goodIndexList(10)]; 

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
    for cond=1:10
        clear oriChangeDegThisCondition selectedPos
        oriChangeDegThisCondition=orientationChangeDeg(goodList{cond});
        selectedPos=find(oriChangeDegThisCondition==uniqueOrientationChangeDeg(ori));
        if isempty(selectedPos)
            firingRate{ori}{cond}=[];
            PSDData{ori}{cond}=[];
            alphaData{ori}{cond}=[];
            gammaData{ori}{cond}=[];
            disp([attCueList{cond} ' Stim Repeats ' num2str(0)]);
            continue
        end
        analyzedLFPData{ori}{cond}=lfpData{cond}.segmentedLFPData(eList,selectedPos,pos);
        analyzedSpikeData{ori}{cond}=spikeData{cond}.segmentedSpikeData(eList,selectedPos);
        disp([attCueList{cond} ' Stim Repeats ' num2str(size(analyzedLFPData{ori}{cond},2))]);
        for k=1:length(eList)
            firingRate{ori}{cond}(k,:)=getSpikeCounts(analyzedSpikeData{ori}{cond}(k,:),timeRange)/diff(timeRange);
            [PSDData{ori}{cond}(k,:,:),freqValsMT]= mtspectrumc(squeeze(analyzedLFPData{ori}{cond}(k,:,:))',params);
        end
        if ~isequal(unique(diff(freqVals)),unique(diff(freqValsMT)))
            error('Frequency values does not match')
        end
        alphaData{ori}{cond}=squeeze(sum(PSDData{ori}{cond}(:,alphaPos,:),2));
        gammaData{ori}{cond}=squeeze(sum(PSDData{ori}{cond}(:,gammaPos,:),2));
    end
end
save(fileToSave,'firingRate','PSDData','alphaData','gammaData','RightArrayPos','LeftArrayPos','attCueList','uniqueOrientationChangeDeg');
end



function [eList,RightArrayPos,LeftArrayPos]=getGoodElectrodes(fileNameString,folderSourceString)

suaCutoff = 3;
% Get sorting rating
sortRatings=sortRating(fileNameString,folderSourceString);
sua=intersect(find(sortRatings>0),find(sortRatings<=suaCutoff));
mua=find(sortRatings>suaCutoff);
eList=union(sua,mua,'sorted');
[~,~,electrodeArrayPos]=electrodePositionOnGridMayo(1,fileNameString);
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