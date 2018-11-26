function [FFC,SFC,PCFF,AmpCorr,alphaCorr,ssvepCorr,Rsc]=getCoherencyAndCorrelation(lfpData,spikeData,electrodePairs,timeRange,alphaPos,ssvepPos)

 timeVals=lfpData{1}.timeVals;
 Fs = round(1/(timeVals(2)-timeVals(1)));
 pos = find(timeVals>=timeRange(1),1)+ (0:diff(timeRange)*Fs-1);
 
 

% Parameters for coherency analysis
params.tapers   = [1 1];
params.pad      = -1;
params.Fs       = 2000;
params.fpass    = [0 200];
params.trialave = 0; 
FFC=cell(1,4); SFC=cell(1,4); PCFF=cell(1,4); AmpCorr=cell(1,4); Rsc=cell(1,4); 
 for array=1:3
    clear FFCTmp SFCTmp PCFFTmp AmpCorrTmp alphaCorrTmp ssvepCorrTmp RscTmp
    for i=1:5
        fftAmp=abs(fft(lfpData{i}.segmentedLFPData(:,:,pos),[],3));
        for j=1:size(electrodePairs{array},1)
            lfp1= squeeze(lfpData{i}.segmentedLFPData(electrodePairs{array}(j,1),:,pos))';
            lfp2= squeeze(lfpData{i}.segmentedLFPData(electrodePairs{array}(j,2),:,pos))';
            spk1= spikeData{i}.segmentedSpikeData(electrodePairs{array}(j,1),:);
            spk2= spikeData{i}.segmentedSpikeData(electrodePairs{array}(j,2),:);
            spkBin=convertSpikeTimes2Bins(spk2,timeRange,1000/Fs);
            spkCnt1=getspkcount(spk1,timeRange);
            spkCnt2=getspkcount(spk2,timeRange);
            amp1=squeeze(fftAmp(electrodePairs{array}(j,1),:,:));
            amp2=squeeze(fftAmp(electrodePairs{array}(j,2),:,:));
            [~,~,sff12,sff1,sff2,~]=coherencyc(lfp1,lfp2,params);
            [~,~,ssf12,ssf1,ssf2,~]=coherencycpb(lfp1,spkBin,params);
            
            FFCTmp(i,j,:)= abs(mean(sff12,2)./sqrt(mean(sff1,2).*mean(sff2,2))); %#ok<AGROW>
            SFCTmp(i,j,:)= abs(mean(ssf12,2)./sqrt(mean(ssf1,2).*mean(ssf2,2)));
            
            coherenceFF= sff12./sqrt(sff1.*sff2);
            PCFFTmp(i,j,:)= abs(mean(coherenceFF,2)); %#ok<AGROW>
            
            RscTmp(i,j)=corr(spkCnt1',spkCnt2');
            
            AmpCorrTmp(i,j,:)=diag(corr(amp1,amp2));
            alphaCorrTmp(i,j)=mean(AmpCorrTmp(i,j,alphaPos),3);
            ssvepCorrTmp(i,j)=AmpCorrTmp(i,j,ssvepPos);
           
        end
    end
    if array==3
        array=array+1;
    end
    FFC{array}=FFCTmp;
    SFC{array}=SFCTmp;
    PCFF{array}=PCFFTmp;
    AmpCorr{array}=AmpCorrTmp;
    alphaCorr{array}=alphaCorrTmp;
    ssvepCorr{array}=ssvepCorrTmp;
    Rsc{array}=RscTmp;
  end    
end

            
function spkcount=getspkcount(X,timeRange)

for i=1:length(X)
    spkcount(i)=length(intersect(find(X{i}>=timeRange(1)),find(X{i}<timeRange(2))));
end

end
            
            