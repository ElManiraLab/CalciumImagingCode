function [outwave] = roi_TSave(Y,roimask)
% [outwave] = roi_TSaverage(Y,roimask). Extracts timeserie of the input
% movie (excluding last frame, in case that the background was included)

% INPUT: 'Y' 3D data matrix (movie: x*y*frames); 'roimask', 2D logic;

% OUTPUT: 'wave' of the ROI timeserie, that is the average value of all the pixels
% in the ROI (roimask)

Nframes= size(Y,3)-1; % substract last frame
outwave = zeros(1, Nframes); % vector of length the nº of frames in data (Y) 

for frame=1:Nframes

        ROIvalues = Y(:,:,frame).*roimask;
        ROIvalues = fillmissing(ROIvalues,'constant', 0); 
        outwave(frame) = sum(ROIvalues(:))/sum(roimask(:)); %sum(ROImask) gives the number of pixels included in the ROI
end

end

%%  Created: 06/12/21 (from old code)
% old function, last debugged on: 21/10/20

% Updated: 07/02/21
