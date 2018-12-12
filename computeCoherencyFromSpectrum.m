% used to calculate the PLF ,ppc values using the spectrum values obtained
% using coherencycpb/coherencyc function;modified from code used in 1/f
% project. 

function [ppc,PLF,without_phaseCovariation,traditional]= computeCoherencyFromSpectrum(s1,s2,s12,params)

% Traditional Coherency
if params.trialave == 0
ms12 = squeeze(mean(s12,3)); ms1 = squeeze(mean(s1,3)); ms2 = squeeze(mean(s2,3));
traditional=abs(ms12./sqrt(ms1.*ms2));


% Coherency without amplitude covariations (PLF)-amplitude normalized to 1;
if params.trialave==1 % for SFC
    product = s12./sqrt(s1.*s2);
    PLF = abs(product);
else
    product = s12./sqrt(s1.*s2);
%     product(isnan(product))=0;
%      PLF = abs(nanmean(product,2));
     PLF = abs(mean(product,2));
end 
% meanPhaseDiff = circ_mean(angle(product),[],2);
% stdPhaseDiff  = circ_std(angle(product),[],[],2);

% Coherency without phase covariations
ms12 = mean(abs(s12),2); ms1 = mean(s1,2); ms2 = mean(s2,2);
without_phaseCovariation= ms12./sqrt(ms1.*ms2);

% PPC
input = (s12./abs(s12))';
siz = size(input);
n = siz(1);
outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
ppc  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way

end