% Save data for Hits versus Miss analysis.

oriChangeList = [2 3];
tpStr = '_TargetOnset'; timePeriod = [-0.5 0];
tapers = [2 3];
alphaRangeHz = [8 12]; gammaRangeHz = [42 78]; lineNoiseFreqHz= 60;
        
% folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
folderSourceString = 'E:\Mayo';

folderNameSave = fullfile(folderSourceString,'Data','savedDataSummary');
makeDirectory(folderNameSave);

fileNameStr = 'dataOri_';
for i=1:length(oriChangeList)
    fileNameStr = cat(2,fileNameStr,num2str(oriChangeList(i)));
end

fileNameStr = cat(2,fileNameStr,tpStr,num2str(timePeriod(1)),'_',num2str(timePeriod(2)),...
    '_tapers', num2str(tapers(1)),num2str(tapers(2)),'_alpha',num2str(alphaRangeHz(1)),'_',num2str(alphaRangeHz(2)),...
    '_gamma',num2str(gammaRangeHz(1)),'_',num2str(gammaRangeHz(2)));
fileNameSave = fullfile(folderNameSave,[fileNameStr '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neutralTrialFlag = 1; 
typeList = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'H0N'} {'H1N'} {'M0N'} {'M1N'}];
[fileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

count=1;
numConditions=length(typeList);
numOrientations=6;
numSessions=25;     % pacu 54 is discarded
numElectrodes=96;

spikeData = cell(numSessions,numConditions,numElectrodes);
alphaData = cell(numSessions,numConditions,numElectrodes);
gammaData = cell(numSessions,numConditions,numElectrodes);

for m=1:length(monkeyNameList)
    for i=1:length(fileNameStringList{m})
        
        disp([count m i]);
        fileNameString = fileNameStringList{m}{i};
        folderName=fullfile(folderSourceString,'Data','segmentedData',fileNameString);

        [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString,neutralTrialFlag);

        for k=1:numConditions
            
            %%%%%%%%%%%%%%%%%%%% Get good positions %%%%%%%%%%%%%%%%%%%%%%%
            x = orientationChangeDeg(goodIndexList{k});
            goodPos = [];
            for j=1:length(oriChangeList)
                goodPos = cat(2,goodPos,find(x==uniqueOrientationChangeDeg(oriChangeList(j))));
            end
            goodPos = sort(goodPos);
            
            % LFP
            lfpDataTMP= load(fullfile(folderName,[fileNameString typeList{k} tpStr '_LFP']));
            timeVals = lfpDataTMP.timeVals;
            timePos = intersect(find(timeVals>=timePeriod(1)),find(timeVals<timePeriod(2)));
            
            % MT analysis
            Fs              = round(1/(timeVals(2)-timeVals(1)));
            params.tapers   = tapers;
            params.pad      = -1;
            params.Fs       = Fs;
            params.fpass    = [0 100];
            params.trialave = 0;
                        
            % Spike
            spikeDataTMP = load(fullfile(folderName,[fileNameString typeList{k} tpStr '_Spikes']));
            
            for j=1:numElectrodes
                
                % LFP
                [psdTMP,freqVals] = mtspectrumc(squeeze(lfpDataTMP.segmentedLFPData(j,goodPos,timePos))',params);
                alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
                lineNoisePos = find(freqVals==lineNoiseFreqHz);
                gammaPos = setdiff(intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<=gammaRangeHz(2))),lineNoisePos);

                alphaData{count,k,j} = sum(psdTMP(alphaPos,:),1);
                gammaData{count,k,j} = sum(psdTMP(gammaPos,:),1);
                
                % Spike
                spikeData{count,k,j} = getSpikeCounts(spikeDataTMP.segmentedSpikeData(j,goodPos),timePeriod);
            end
        end
        count=count+1;
    end
end

save(fileNameSave,'spikeData','alphaData','gammaData');