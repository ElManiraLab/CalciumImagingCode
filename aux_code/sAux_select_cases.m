condit = [1 2]; %@ SET
iti_thresh = 45; %@ SET seconds (minimum iti)

% run this line to select the trials with the condition, meeting the iti
% requirement
idx = intersect(choose_condidx(VSDI.condition, condit), find(VSDI.iti>iti_thresh)); 

