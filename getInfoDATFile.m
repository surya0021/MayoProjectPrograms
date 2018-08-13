function [DAT2,CueOnCode,targetOnTimeMS,orientationChangeDeg,saccadeTimeMS] = getInfoDATFile(DAT)
i=1; continueFlag=1;
while(continueFlag)
    tnum = (['t', num2str(i)]);
    if isfield(DAT,tnum)
        DAT2{i} = DAT.(tnum); %#ok<*AGROW>
        if isfield(DAT2{i},'cueOn')
            CueOnCode(i) = 1;
        else
            CueOnCode(i) = 0;
        end
        
        if isfield(DAT2{i}.trial.data,'targetOnTimeMS')
            targetOnTimeMS(i) = DAT2{i}.trial.data.targetOnTimeMS; % From stimulus onset
            orientationChangeDeg(i) = DAT2{i}.trial.data.change;
        else
           targetOnTimeMS(i) = 0;
           orientationChangeDeg(i) = 0;
        end
        
        if isfield(DAT2{i},'saccade')
            saccadeTimeMS(i) = DAT2{i}.saccade.timeMS - DAT2{i}.trialStart.timeMS;
        else
           saccadeTimeMS(i) = 0;
        end
        
        i=i+1;
    else
        continueFlag=0;
    end
end
end