function [cropped] = roi_crop(movie3D, cropmask)
%INPUT:
% 'movie3D': dimensions: (x*y*frame)
% 'cropmask': 2D logic array (mask) to apply at each movies frame

% APPLY CROPMASK TO EACH MOVIE FRAME
nframe = size(movie3D, 3);
cropped = zeros(size(movie3D));
crop3D = repmat(cropmask,  1,1,nframe);
cropped = movie3D .* crop3D; 

end

%% Created: 06/02/21 from older code
% Updated: 