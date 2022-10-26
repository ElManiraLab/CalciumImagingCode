function [output2] = filter_cubicBV(input)

% SPATIAL + TEMPORAL FILTER that mimics the Cubic BrainVision filter. 
% INPUT : 3D movie with background (only the spatial component will be
% applied)
% OUTPUT: filtered 3D movie

% NOTE on the temporal filtering: this is the one used by the program but it could be
% preferable to use a filter that smoothes in both directions

% -----------------
% Brainvision definition: "Cubic filter processing" performs a "3 x 3," 2D filter processing and temporal filter processing. 
% The value defined as D(t, x, y) at the coordinates of (x, y) on the frame number "t" will be is calculated with following formula.

% V(t,x,y) = D(t,x,y)/2 + D(t,x-1,y)/8 + D(t,x+1,y)/8 + D(t,x,y-1)/8 + D(t,x,y+1)/8
% D(t,x,y) = V(t-1,x,y)/4 + V(t,x,y)/2 + V(t+1,x,y)/4

n_frame = size(input,3); 

% TEMPORAL SMOOTHING

    for t= 2:n_frame-2 %we exclude extreme frames (and exclude an extra one to leave the background out of the comutations)
    output(:,:,t) = input(:,:,t)/2 +  input(:,:,t-1)/4 + input(:,:,t+1)/4;  
    end

    % extreme frames are smoothed in just one direction
    output(:,:,1) = input(:,:,1)/2 + input(:,:,2)/2;
    output(:,:,n_frame) = input(:,:,end-1)/2 + input(:,:,end)/2; %'end-1' es el último frame de actividad

%SPATIAL SMOOTHING (according to formula)

    kernel = [0 1/8 0; 1/8 1/2 1/8; 0 1/8 0];

    for fr = 1:n_frame
    frame = squeeze(input(:,:,fr));
    kernel = kernel / sum(kernel(:));
    localAverageImage = conv2(frame, kernel, 'same');
    output2(:,:,fr) = localAverageImage; 
    end

end

%% Created: 17/02/21 (from Gent2)