% saveGoodElectrodes
populationType = 'All';

folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Surya_MayoProject';
folderNameSave = fullfile(folderSourceString,'Data','savedDataSummary');
makeDirectory(folderNameSave);

[fileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

count=1;

for m=1:length(monkeyNameList)
    for i=1:length(fileNameStringList{m})
        
        disp([count m i]);
        fileNameString = fileNameStringList{m}{i};
        
        electrodeArrayList{count} = getGoodElectrodes(fileNameString,folderSourceString,populationType);
        count=count+1;
    end
end

fileNameSave = fullfile(folderNameSave,['electrodeArrayList' populationType '.mat']);
save(fileNameSave,'electrodeArrayList');