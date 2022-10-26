
function [movie_thres, alphachan, thresh_used] = movie_noisethresh(movie,baseline,SDfactor, return_thresh)

% Set to zero all the values below the 'noise
% threshold'. The 'noise threshold' is computed at each pixel, and is
% defined as SDfactor times the STD of the signal in the baseline.

% INPUT:
% 'movie': 3Dmovie
% baseline: array of idx of frames to be considered reference baseline
% (noise)
% SDfactor: single value
% 'return_thresh' 0 or 1

% OUTPUT:
% 'movie_thres' - thresholded movie
% 'alphachan' - alpha channel to apply transparency

movie_thres=zeros(size(movie)); % preallocate

for rowi= 1:size(movie,1)
    for coli= 1:size(movie,2)
        s  = squeeze(movie(rowi, coli,:));
        noiseSD = std(s(baseline));
        noisemean = mean(s(baseline));
        
        s=s-noisemean; % y-shifting the signal allows to threshold directly with the SD 
        zero_out = find(abs(s)<=SDfactor*noiseSD); 
        s_th= s;
        s_th(zero_out) = 0; % points where |s| is less than SDfactor*SD are "sub-threshold" and we set them to zero.

        movie_thres(rowi,coli,:) = s_th; 
        % the last input to SDthresh is how many times away from the standard deviation the signal must be before we say the signal isn't noise
        
        if return_thresh
            thresh_used(rowi,coli) = SDfactor*noiseSD; % threshold used for each timepoint
        else
            thresh_used= [];
        end
        
    end
    
    alphachan = movie_thres==0 ;
    
end

%% Created: 07/02/21 (copied from Gent2's code)
