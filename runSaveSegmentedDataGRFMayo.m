% This is the main program for doing all data extraction of GRF data.

% fileNameString: string indicating the fileName.
% We assume that the following files are available at
% folderSourceString\Data\extractedData
% 1. fileNameString_DAT
% 2. fileNameString_extractedTrialsLFP
% 3. optionally, fileNameString_Codes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileNameString = 'pwmp29A001';
folderSourceString='F:\Projects\MayoProject\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 - Get Lablib LL Data
saveLLDataGRFMayo(fileNameString,folderSourceString);

% Step 2 - extract stimulus information in a useful format
extractStimInfoGRFMayo(fileNameString,folderSourceString,0);
extractStimInfoGRFMayo(fileNameString,folderSourceString,1);

% Step 3 - generate 'parameterCombinations' that allows us to find the useful combinations
getDisplayCombinationsGRFMayo(fileNameString,folderSourceString,0); % side0
getDisplayCombinationsGRFMayo(fileNameString,folderSourceString,1); % side1

% Step 4 - segment and save data
saveSegmentedDataGRFMayo(fileNameString,folderSourceString);