function [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='C:\Supratim\Projects\MayoProject\';       end

folderNameIn = fullfile(folderSourceString,'Data','extractedData');

fileNameCDS = [fileNameString '_Codes'];
CDS = load(fullfile(folderNameIn,fileNameCDS));
CDS = CDS.(fileNameCDS);

fileNameDAT = strcat(fileNameString, '_DAT');
DAT = load(fullfile(folderNameIn,fileNameDAT));
DAT = DAT.(fileNameDAT);

goodIndexList = getGoodIndices(CDS,DAT);
[~,~,~,orientationChangeDeg,~] = getInfoDATFile(DAT);

uniqueOrientationChangeDeg = unique(orientationChangeDeg);

numConditions = length(goodIndexList);
numOrientations = length(uniqueOrientationChangeDeg);

responseMatrix = zeros(numConditions,numOrientations);
for i=1:numConditions
    x = orientationChangeDeg(goodIndexList{i});
    for j=1:numOrientations
        responseMatrix(i,j) = length(find(x==uniqueOrientationChangeDeg(j)));
    end
end

perCorrect = zeros(5,numOrientations); % order: 0V, 1V, 0I, 1I, N
for i=1:4
    perCorrect(i,:) = responseMatrix(i,:) ./ (responseMatrix(i,:)+responseMatrix(i+4,:));
end
perCorrect(5,:) = responseMatrix(9,:) ./ (responseMatrix(9,:)+responseMatrix(10,:));
end