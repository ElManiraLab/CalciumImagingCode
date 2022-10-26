function [registmov] = register_wrap(rawmov, ref_frame,selected_trials)
% 

% INPUT
% 'rawmov' - unregistered 4D movies with dimensions: x*y*time*trials
% 'ref_frame' - 2Dmatrix containing reference frame for registering
% 'selected_trials' -  idx of trials that contain movie (non-NaN VSDI.condition) ?
% that are not empty

%OUTPUT
% 'registmov' - movie aligned respect to the ref_frame

% INPUT CONTROL 
if nargin  <3
selected_trials = 1:size(rawmov,4);
end

% Check that selected trials is a row (for the loop)
if size(selected_trials,1) > size(selected_trials,2)
    selected_trials = selected_trials';
end

% END OF INPUT CONTROL

%% 'imregister' CONFIGURATION

% Registration configuration: MONOMODAL
    % Get a configuration suitable for registering images from different
    % sensors.
    [optimizer, metric] = imregconfig('monomodal'); 


% OLD CODE:  mutimodal, specifying parameters ----
%     % Get a configuration suitable for registering images from different
%     % sensors.
%     [optimizer, ~] = imregconfig('multimodal'); 
%  
%     % Tune the properties of the optimizer to get the problem to converge
%     % on a global maxima and to allow for more iterations.
%     optimizer.InitialRadius = 0.009;
%     optimizer.Epsilon = 1.5e-4;
%     optimizer.GrowthFactor = 1.01;
%     optimizer.MaximumIterations = 500;
% 
% metric = registration.metric.MattesMutualInformation;
% -------

%% ALIGN ALL MOVIES, MATCHING THE FIRST FRAME TO THE REFERENCE FRAME (ref_frame)
% import all movies and align respect to the reference frame
for idx =  selected_trials
    %Get the transformations that matches current trial's background with the reference background
    [tform]= imregtform(squeeze(rawmov(:,:,1,idx)),ref_frame, 'rigid', optimizer, metric); 
    Rfixed = imref2d(size(ref_frame));% coordinateds

    % Apply the transformations to the rest of the frames
    for fri = 1:size(rawmov,3)
        
%    registdata(:,:,fri)= imregister_mod(rawdata(:,:,fri),frame_ref, 'rigid', optimizer, metric); 
     registmov(:,:,fri,idx) = imwarp(squeeze(rawmov(:,:,fri,idx)),tform,'OutputView',Rfixed, 'SmoothEdges', false);
    end
    
    display(idx)
end
end