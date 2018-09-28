% This function gives the detection performance(%correct) in one session of
% left,right or both the visual fields
% side: 0-left,1-right,2-Both.


function percorrect=getBehavior(fileNameString,side)
cue=[{'V'} {'I'}];
for attCond=1:3
    if attCond==3
        conditionString=[{'HN'} {'MN'}];
    else
        if side==0
            conditionString=[{['H0' cue{attCond}]} {['M0' cue{attCond}]}];
        elseif side==1
            conditionString=[{['H1' cue{attCond}]} {['M1' cue{attCond}]}];
        elseif side==2
            conditionString=[{['H0' cue{attCond}]} {['M0' cue{attCond}]} {['H1' cue{attCond}]} {['M1' cue{attCond}]}];
        end
    end
    hit=zeros(1,6);
    miss=zeros(1,6);
    for j=1:length(conditionString)
        for k=1:6
            goodPos{k}=getGoodTrials(fileNameString,conditionString{j},k); %#ok<AGROW>
            nPos(k)=length(goodPos{k}); %#ok<AGROW>
        end
        if bitget(j,1)
            hit=hit+nPos;
        else
            miss=miss+nPos;
        end
    end
    percorrect(attCond,:)=(hit./(hit+miss))*100; %#ok<AGROW>
end
end