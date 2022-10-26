function [sharks] =  find_sharks (GS, VSDI, SDthresh, plot_on)

% Find trials with outliers (mostly shark-like) responses according to the
% SDthreshold input. It performs it both in baseline (using all trials) and
% post-stimulus (condition-specific-wise)


% INPUT
%  'GS': 2d array (timepoints, trials). Contains the whole-brain
%   pixels-average (GS) for each trial
% 'VSDI' structure


% OUTPUT
% 'trial_idx' - idx of trials to reject


% input control
if nargin <4
    plot_on = 0;
end
% end of input control

ptp = dsp.PeakToPeak();


[~, idxOnset] = min(abs(VSDI.times - VSDI.design.Sonset)); 
idxend = length(VSDI.times); 

for triali = VSDI.trials_in'
sm_GS = medfilt1(GS(:,triali), 15); sm_GS =filter_Tcnst(sm_GS,10);
maxpeakPRE(triali)= ptp(sm_GS(1:idxOnset));
maxpeakPOST(triali)= ptp(sm_GS(idxOnset:idxend));
end

peak_meanPRE = mean(maxpeakPRE(VSDI.trials_in));
peak_stdPRE = std(maxpeakPRE(VSDI.trials_in));

peak_mean0 = mean(maxpeakPOST(intersect(find(VSDI.condition ==0), VSDI.trials_in)));
peak_std0 = std(maxpeakPOST(intersect(find(VSDI.condition ==0), VSDI.trials_in)));
peak_mean1 = mean(maxpeakPOST(intersect(find(VSDI.condition ==1), VSDI.trials_in)));
peak_std1 = std(maxpeakPOST(intersect(find(VSDI.condition ==1), VSDI.trials_in)));
peak_mean2 = mean(maxpeakPOST(intersect(find(VSDI.condition ==2), VSDI.trials_in)));
peak_std2 = std(maxpeakPOST(intersect(find(VSDI.condition ==2), VSDI.trials_in)));

sharks = [];
SDthresh = 2;
% absthresh = 0.3; % absolute threshold the signal cannot pass

maxpeakPRE = maxpeakPRE - peak_meanPRE;
for triali = VSDI.trials_in'
if maxpeakPRE(triali) > peak_std1*SDthresh  
    sharks = [sharks triali] ; 
end
end

% condition 0
maxpeak0 = maxpeakPOST - peak_mean0;
%maxpeakPOST includes all-condition trials, but then just Condition-specific are considered in the loop
for triali = intersect(find(VSDI.condition ==0), VSDI.trials_in)'

% sm_GS = medfilt1(GS(:,triali), 15); sm_GS =filter_Tcnst(sm_GS,10);
% crossidx = sm_GS;
% maxpeak 
if maxpeak0(triali) > peak_std0*SDthresh 
    sharks = [sharks triali] ; 
end
end

% condition 1
maxpeak1 = maxpeakPOST - peak_mean1; 
for triali = intersect(find(VSDI.condition ==1), VSDI.trials_in)'
if maxpeak1(triali) > peak_std1*SDthresh  
    sharks = [sharks triali] ; 
end
end

% condition 2
maxpeak2 = maxpeakPOST - peak_mean2;
for triali = intersect(find(VSDI.condition ==2), VSDI.trials_in)'
if maxpeak2(triali) > peak_std2*SDthresh 
    sharks = [sharks triali] ;
end
end

sharks = sort(sharks)'; 

if plot_on
figure
plot(GS(:,sharks), 'color','#97978F'); hold on;
plot(GS(:,setdiff(VSDI.trials_in,sharks)), 'color','#FBC4BF'); 
end
end

