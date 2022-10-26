function [outwave] = roi_TSaverage(Y,roimask)
% [outwave] = roi_TSaverage(Y,roimask). Extracts timeserie 

% INPUT: 'Y' 3D data matrix (movie: x*y*frames); 'roimask', 2D logic;

% OUTPUT: 'wave' of the ROI timeserie, that is the average value of all the pixels
% in the ROI (roimask)

Nframes= size(Y,3);
outwave = zeros(1, Nframes); % vector of length the nº of frames in data (Y) 

for frame=1:Nframes

        ROIvalues = Y(:,:,frame).*roimask;
        outwave(frame) = sum(ROIvalues(:))/sum(roimask(:)); %sum(ROImask) gives the number of pixels included in the ROI
end

end

%%  Created: 06/12/21 (from old code)
% old function, last debugged on: 21/10/20
