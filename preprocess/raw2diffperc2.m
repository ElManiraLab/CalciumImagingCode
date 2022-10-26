function [diffpdata] = raw2diffperc2(raw3Dmatrix, baseframes)
%RAW (absolute intensity) MOVIES TO DIFFERENTIAL PERCENT VALUE
%-[100*(F-F0)/F0] MOVIES, BEING F0 THE MEAN OF THE BASERANGE FRAMES (idx)

% For each pixel, it calculates differential values respect to the F0

% INPUT   'raw3Dmatrix' raw data MOVIE matrix (hpix, wpix, n_frames) 
%         'baseframes': idx of baseframe to calculate F0.

% OUTPUT   'diffpdata'  diffenrential values MOVIE(=3D)  matrix -100*(F-F0)/F0; raw-values
% background (frame0 from the rawmovie) is stored in an extra frame in diffpdata(:,:,n_frames+1)               

% example VSDmov.data(:,:,:,1) = raw2diffperc2(raw3Dmatrix, 1);

imdim1 = size(raw3Dmatrix,1); %1st dimension of the movie frame
imdim2 = size(raw3Dmatrix,2); 
n_frames = size(raw3Dmatrix,3); % nframes

% Compute baseline (F0) mean value:
selframes = raw3Dmatrix(:,:,baseframes);
F0 = mean(selframes,3);

% Replacing 0 for a small value to avoid inf or -inf
F0(F0==0)=1.5259e-05;

% Compute 100*(F-F0)/F0 for each value
diffpdata = zeros(imdim1,imdim2,n_frames+1);

for framei= 1:n_frames
    diffpdata(:,:,framei)= ( 100*(raw3Dmatrix(:,:,framei)-F0) )./F0;       
end

% Invert Sign
diffpdata = -diffpdata;

% Store background in the last frame
diffpdata(:,:,n_frames+1) = raw3Dmatrix(:,:,1);

% Don't remember why the next line of code; let's try without it (there
% should be no NaN's values)
% diffpdata = fillmissing(diffpdata,'constant',0);
end 

%% Created: 03/02/2021 (from raw2diffperc)
% Updated: 