function [output] = filter_median(input, pix)
% [output] = filter_median(input)
% pix*pix MEDIAN FILTER OF THE IMAGE (2D input) OR OF THE FRAMES OF A
% MOVIE (3Dinput). If 'pix' is not input, 3*3 window is used

% INPUT: image or movie. If it is a movie, the median filter is applied in
% the first two dimensions
% OUTPUT: filtered image or movie

% Performs "pix x pix" pixels average filter of 2D filter processing.

% If The value defined as D(t, x, y) on the coordinates of (x, y) at the frame number "t" will be replaced with the averaged value - in the area "P x P".
% If input is an image (dim=2), returns an image filtered.
% If the input is a movie (3D), ut returns a movie in which each frame has
% been filtered independently

% Input control
if nargin <2
pix = 3; 
end
% end of input control

%% 
data_dim = length(size(input));

if data_dim ==2 
output = medfilt2(input, [pix pix]); %performs 3x3pixels median filter
end

if data_dim ==3
n_frame = size(input,3); %nº time points

for fr = 1:n_frame % apply the filter to each frame
frame = squeeze(input(:,:,fr));
localAverageImage = medfilt2(frame, [pix pix]);
output(:,:,fr) = localAverageImage; 
end

end

end

%% Created: 13/02/21 (from Gent's code)
% Updated:

% DAVID: lo que que se pueda ajustar la ventana lo he puesto nuevo;
% a ver si puedes comprobarlo