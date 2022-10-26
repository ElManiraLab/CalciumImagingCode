
% >>>> Manually set the condition on the last column: 0 = no_tone; 1 =
% 1000hz tone; 2 = 2000Hz tone; NaN = not classifyed

% Fill which rows correspond to the different conditions...
range.tone1000 = [11:20,41:50,70:79,100:109,130:139];
range.no_tone = [21:30,51:59, 80:89,110:119,140:149];
range.tone2000 = [31:40,60:69,90:99, 120:129,150:159];

% ... and run this code
for ii = 1:length(VSDI.list)
VSDI.list(ii).condition = NaN;
end

for ii = range.no_tone
VSDI.list(ii).condition = 0;
end

for ii = range.tone1000
VSDI.list(ii).condition = 1;
end

for ii = range.tone2000
VSDI.list(ii).condition = 2;
end

% extract the trials indexes (condition ~=NaN)
for ii = 1:length(VSDI.list)
condition(ii) = VSDI.list(ii).condition;
ref(ii)= str2num(VSDI.list(ii).name(1:end-5));
end
condition = condition'; ref = ref';

VSDI.condition = condition;
VSDI.trialref = ref;
clear ref condition range filelist* labels
% % % 
