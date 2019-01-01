% used to calculate the PLF ,ppc values using the spectrum values obtained
% using coherencycpb/coherencyc function; modified from code used in 1/f
% project

function [ppc,PLF,without_phaseCovariation,traditional]= computeCoherencyFromSpectrum(s1,s2,s12,params)

% Traditional Coherency (Calculate Coherency on trial averaged data)
if params.trialave == 0
    ms12 = squeeze(mean(s12,2)); ms1 = squeeze(mean(s1,2)); ms2 = squeeze(mean(s2,2)); 
    traditional = abs(ms12./sqrt(ms1.*ms2)); % Coherency formula
elseif params.trialave == 1 
    traditional = abs(S12./sqrt(S1.*S2)); % Coherency formula 
end

% Coherency without amplitude covariations (PLF)-amplitude normalized to 1;
% (Calculate Coherency on each trial and take trial average later)- Phase
% Coherence
if params.trialave==0 
     product = s12./sqrt(s1.*s2);
    %      product(isnan(product))=0;
    %      PLF = abs(mean(product,2));
     PLF = abs(nanmean(product,2));
    % meanPhaseDiff = circ_mean(angle(product),[],2);
    % stdPhaseDiff  = circ_std(angle(product),[],[],2);
elseif params.trialave==1
    product = s12./sqrt(s1.*s2);
    PLF = abs(product);
end 

% Coherency without phase covariations- Amplitude Coherence
if params.trialave==0 
    ms12 = mean(abs(s12),2); ms1 = mean(s1,2); ms2 = mean(s2,2);
    without_phaseCovariation= ms12./sqrt(ms1.*ms2);
elseif params.trialave==1
    without_phaseCovariation = abs(s12)./sqrt(s1.*s2);
end

% PPC works only for non-trial averaged dataset because it is a trialwise
% measure
% PPC- adapted from fieldtrip/connectivity/ft_connectivity_ppc.m computes 
% pairwise phase consistency  from a data-matrix containing a cross-spectral 
% density. This implements the method described in Vinck M, van Wingerden M, 
% Womelsdorf T, Fries P, Pennartz CM.
% The pairwise phase consistency: a bias-free measure of rhythmic neuronal
% synchronization. Vinck et al. Neuroimage. 2010

input = (s12./abs(s12))'; % normalize the cross-spectrum
siz = size(input);
n = siz(1);
if n>1
    outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
    ppc  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
else
    error('computation of PPC requires >1 trial, please feed all trial dataset into computeCoherencyFromSpectrum program')
end
end