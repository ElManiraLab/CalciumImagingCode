function [raw4D, times] =  assemble4Draw(VSDI,folder,stime)
% takes included individual movies from the 'folder' and assembles a

% INPUT
% VSDI structure with:
%   'VSDI.condition' filed (sets the inclusion criterion: non-NaN value in VSDI.condition)
%   'VSDI.trialref': to get the names of the .dml files to import
%   'stime': to calculate the time base

% OUTPUT
% movies4D (x*y*time*trials)
% In time dimension: time=1 >> frame0(0ms=background);
% time=2:681 >> frames1-to-end(stime-stime*nframes ms);

for triali = 1:length(VSDI.condition)
    if ~isnan(VSDI.condition(triali,1)) %Check whether the first condition is not NaN;
        % if the first condition column is not Nan, it's enough to know
        % that it's a valid trial
        load(fullfile(folder,[num2str(VSDI.trialref(triali)) 'A.mat']))
        data_act = permute(data_act, [2 3 1]);
        data_act = double(data_act);
        raw4D(:,:,:,triali) =  data_act;
    end
end

% Positive values only:
minval = min(raw4D(:));
raw4D = raw4D +abs(minval); 

% get timebase
nframes = size(raw4D,3);
times = (0:nframes-1)*stime;
times = times';

end

%% Created: 01/02/2021
% Updated: 03/02/2021