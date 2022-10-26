function [output] =CALCIUMimg(action, CALCIUM, CALCIUM_feature,list,nfish)
% PERFORMS BASIC LOAD/SAVE FUNCTIONS REFERENCE TOCALCIUMimg
% action = 'save' or 'load' or 'savemovie' or 'loadmovie'

% Input and output will depend on the case ('action')
% Use according to 'action':
%   [~]=CALCIUMimg('save', CALCIUM) or CALCIUMimg('save', CALCIUM)-uses internal CALCIUM.ref to save in appropiate
%	.mat
%   [CALCIUM] =CALCIUMimg('load', nfish)
%  [~]=CALCIUMimg('savemovie', VSDmov, movierefernce) - uses moviereference
%  (~char) to name the matfile

%  [VSDmov]=CALCIUMimg('loadmovie', nfish, moviereference) - uses moviereference

path.Calcium_ourToolbox = pwd ;  % run the code from the inside of the rootpath folder
datapath = fullfile(path.Calcium_ourToolbox(1:end-8),'Calcium_data');


CALCIUMpath = fullfile(datapath,'dataCALCIUM',char(list(nfish,1))); % Saving
moviepath = fullfile(datapath,'datamovies',char(list(nfish,1))); % Saving
wavespath = fullfile(datapath,'datawaves',char(list(nfish,1))); % Saving

expref = 'CALCIUM';

% Input control
switch action
    case 'save'
            if  ~isstruct(CALCIUM) 
            disp('the input is not what expected'); end
    case 'load'
%         assert(mod(CALCIUM, 1) == 0 && , 'input to load must be a single number');
        
        try
            load(fullfile(datapath, 'grouplist.mat'))
        catch 
            warning('fish cannot be load because "grouplist.mat" does not exist')
        end
        
     case 'savemovie'
            if ~exist('CALCIUM_feature') 
                error('input a proper reference name for the movie (as 3rd argument)'); end
end % input control

%% FUNCTION CODE:

switch action
    case 'save'
        CALCIUM = CALCIUM; 
        %saveCALCIUM saves current CALCIUM structure respect to the current rootpath
        pathname = fullfile(CALCIUMpath,[expref '_' CALCIUM.ref '.mat']);
        save(pathname, 'CALCIUM')

    case 'load'
        load(fullfile(datapath, 'grouplist')) %load structure list to take the fish reference
        load(fullfile(CALCIUMpath,[grouplist{CALCIUM},'.mat'])) %load the fish CALCIUM
        disp(strcat (grouplist{nfish}, '_loaded'));
        output= CALCIUM;
        
    case 'savemovie' 
       VSDmov= CALCIUM;
       %saveCALCIUM saves current CALCIUM structure respect to the current rootpath
       pathname = fullfile(moviepath,['CALCIUMmov_',num2str(VSDmov.ref),CALCIUM_feature,'.mat']);
       save(pathname,'VSDmov', '-v7.3')

    case 'loadmovie' 
       load(fullfile(datapath, 'grouplist'))
       fishref = grouplist{CALCIUM}(9:end);
       %saveCALCIUM saves current CALCIUM structure respect to the current rootpath
       movieref = [expref,'Mov_',fishref,CALCIUM_feature,'.mat'];
       load(fullfile(moviepath,movieref))
       output= VSDmov;       
       disp([movieref, '_loaded']);

    case 'savewave'
        VSDroiTS = CALCIUM; 
        %saveCALCIUM saves current CALCIUM structure respect to the current rootpath
        pathname = fullfile(wavespath,[expref 'RoiTS_' num2str(CALCIUM.ref) '.mat']);
        save(pathname, 'VSDroiTS')

    case 'loadwave'
        load(fullfile(datapath, 'grouplist')) %load structure list to take the fish reference
        fishref = grouplist{CALCIUM}(9:end);
        load(fullfile(wavespath,[expref 'RoiTS_',fishref,'.mat'])) %load the fish CALCIUM
        disp(strcat ('ROIs timeseries for fish',grouplist{CALCIUM}, '_loaded'));
        output= VSDroiTS;



        
end %switch
end

% function T = isIntegerValue(X)
% T = (mod(X, 1) == 0);
% end

%% Created: 31/01/2021
% Updated: 08/02/21
