fileNameString = 'paoc04A001';
folderSourceString='F:\Projects\MayoProject\';
folderName = fullfile(folderSourceString,'Data','segmentedData',fileNameString);

x=load(fullfile(folderName,[fileNameString 'H0V_StimOnset'])); % AttLoc 1
y=load(fullfile(folderName,[fileNameString 'H1V_StimOnset'])); % AttLoc 2
locationStr = [{'Right'} {'Left'}];

electrodeArray{1} = [75 73 71 69 67 65 72 70 68 66 79 77 02 81 80 78 76 74]; % Right Array
electrodeArray{2} = [46 50 15 17 09 06 11    90 89 55 56 57 58 52 54 19 13 25]; % Left Array

timeVals=x.timeVals;
blRange = [-0.25 0]; stRange = [0.25 0.5];
Fs = round(1/(timeVals(2)-timeVals(1)));
blPos = find(timeVals>=blRange(1),1)+ (0:diff(blRange)*Fs-1);
stPos = find(timeVals>=stRange(1),1)+ (0:diff(stRange)*Fs-1);

xsBL = 0:1/(diff(blRange)):Fs-1/(diff(blRange));
xsST = 0:1/(diff(stRange)):Fs-1/(diff(stRange));

ssvepPos=6;
for i=1:2
    eList = electrodeArray{i};
    
    subplot(2,3,3*(i-1)+1);
    erp1 = squeeze(mean(mean(x.segmentedLFPData(eList,:,:),1),2));
    erp2 = squeeze(mean(mean(y.segmentedLFPData(eList,:,:),1),2));
    plot(x.timeVals,erp1,'b'); hold on; 
    plot(y.timeVals,erp2,'r');
    legend('Loc0','Loc1');
    xlim([-0.25 0.5]);
    title(['ERP, Array Side: ' locationStr{i}]);
    xlabel('Time from stim onset (s)'); ylabel('Voltage (\muV)');
    
    subplot(2,3,3*(i-1)+2);
    fftst1 = squeeze(mean(log10(abs(fft(x.segmentedLFPData(eList,:,stPos),[],3))),2));
    fftst2 = squeeze(mean(log10(abs(fft(y.segmentedLFPData(eList,:,stPos),[],3))),2));
    fftbl1 = squeeze(mean(log10(abs(fft(x.segmentedLFPData(eList,:,blPos),[],3))),2));
    fftbl2 = squeeze(mean(log10(abs(fft(y.segmentedLFPData(eList,:,blPos),[],3))),2));
    
    dPower1=10*(fftst1-fftbl1); dPower2=10*(fftst2-fftbl2);
    plot(xsST,mean(dPower1,1),'b'); hold on;
    plot(xsST,mean(dPower2,1),'r');
    plot(xsST(ssvepPos),mean(dPower1(:,ssvepPos)),'bo');
    plot(xsST(ssvepPos),mean(dPower2(:,ssvepPos)),'ro');
   
    legend('Loc0','Loc1');
    xlim([0 100]);
    title('Change in FFT amplitude (stim-baseline)');
    xlabel('Frequency (Hz)');
    
    subplot(2,3,3*(i-1)+3);
    dPower1SSVEP = dPower1(:,ssvepPos); dPower2SSVEP = dPower2(:,ssvepPos);
    bar(0,mean(dPower1SSVEP),'b'); hold on;
    bar(1,mean(dPower2SSVEP),'r'); hold on;
    errorbar(0,mean(dPower1SSVEP),std(dPower1SSVEP)/sqrt(length(dPower1SSVEP)),'k');
    errorbar(1,mean(dPower2SSVEP),std(dPower2SSVEP)/sqrt(length(dPower2SSVEP)),'k');
    title('Change in SSVEP amplitude (stim-baseline)');
end