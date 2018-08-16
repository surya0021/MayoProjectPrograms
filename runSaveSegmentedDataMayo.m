[fileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

for m=1:length(monkeyNameList)
    for i=1:length(fileNameStringList{m})
        saveSegmentedDataMayo(fileNameStringList{m}{i});
    end
end