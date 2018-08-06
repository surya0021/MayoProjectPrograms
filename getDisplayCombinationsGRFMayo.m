% This program generates the parameterCombinations variable from stimResults
function parameterCombinations = getDisplayCombinationsGRFMayo(fileNameString,folderSourceString,activeSide,radiusToSigmaRatio)

if ~exist('radiusToSigmaRatio','var');      radiusToSigmaRatio=3;       end

folderNameOut = fullfile(folderSourceString,'Data','segmentedData',fileNameString);
load(fullfile(folderNameOut,['stimResults' num2str(activeSide) '.mat']));
load(fullfile(folderNameOut,['goodStimNums' num2str(activeSide) '.mat']));

% Five parameters are chosen.
% 1. Azimuth
% 2. Elevation
% 3. Sigma, Radius 
% 4. Spatial Frequency
% 5. Orientation

% Parameters index
parameters{1} = 'azimuth';
parameters{2} = 'elevation';
parameters{3} = 'sigma';
parameters{4} = 'spatialFrequency';
parameters{5} = 'orientation'; %#ok<NASGU>

% get Data
aValsAll  = stimResults.azimuth;
eValsAll  = stimResults.elevation;
if isfield(stimResults,'radius')
    sValsAll  = stimResults.radius/radiusToSigmaRatio;
else
    sValsAll  = stimResults.sigma;
end
fValsAll  = stimResults.spatialFrequency;
oValsAll  = stimResults.orientation;

if ~isempty(aValsAll)

    aValsGood = aValsAll(goodStimNums);
    eValsGood = eValsAll(goodStimNums);
    sValsGood = sValsAll(goodStimNums);
    fValsGood = fValsAll(goodStimNums);
    oValsGood = oValsAll(goodStimNums);

    aValsUnique = unique(aValsGood); aLen = length(aValsUnique);
    eValsUnique = unique(eValsGood); eLen = length(eValsUnique);
    sValsUnique = unique(sValsGood); sLen = length(sValsUnique);
    fValsUnique = unique(fValsGood); fLen = length(fValsUnique);
    oValsUnique = unique(oValsGood); oLen = length(oValsUnique);

    % display
    disp(['Number of unique azimuths: ' num2str(aLen)]);
    disp(['Number of unique elevations: ' num2str(eLen)]);
    disp(['Number of unique sigmas: ' num2str(sLen)]);
    disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
    disp(['Number of unique orientations: ' num2str(oLen)]);

    % If more than one value, make another entry with all values
    if (aLen > 1);           aLen=aLen+1;                    end
    if (eLen > 1);           eLen=eLen+1;                    end
    if (sLen > 1);           sLen=sLen+1;                    end
    if (fLen > 1);           fLen=fLen+1;                    end
    if (oLen > 1);           oLen=oLen+1;                    end

    allPos = 1:length(goodStimNums);
    disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen))]);

    for a=1:aLen
        if a==aLen
            aPos = allPos;
        else
            aPos = find(aValsGood == aValsUnique(a));
        end

        for e=1:eLen
            if e==eLen
                ePos = allPos;
            else
                ePos = find(eValsGood == eValsUnique(e));
            end

            for s=1:sLen
                if s==sLen
                    sPos = allPos;
                else
                    sPos = find(sValsGood == sValsUnique(s));
                end

                for f=1:fLen
                    if f==fLen
                        fPos = allPos;
                    else
                        fPos = find(fValsGood == fValsUnique(f));
                    end

                    for o=1:oLen
                        if o==oLen
                            oPos = allPos;
                        else
                            oPos = find(oValsGood == oValsUnique(o));
                        end
                        
                        aePos = intersect(aPos,ePos);
                        aesPos = intersect(aePos,sPos);
                        aesfPos = intersect(aesPos,fPos);
                        aesfoPos = intersect(aesfPos,oPos);
                        
                        parameterCombinations{a,e,s,f,o} = aesfoPos; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % save
    save(fullfile(folderNameOut,['parameterCombinations' num2str(activeSide) '.mat']),'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique');

end
end