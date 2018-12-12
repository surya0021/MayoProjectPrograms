[fileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

folderSourceString='C:\Supratim\Projects\MayoProject\';
instructionTrialFlag=1;

for m=1:length(monkeyNameList)
    for i=1:length(fileNameStringList{m})
        saveSegmentedDataMayo(fileNameStringList{m}{i},folderSourceString,instructionTrialFlag);
    end
end