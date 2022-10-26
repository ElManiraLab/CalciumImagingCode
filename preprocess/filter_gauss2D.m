function [output] = filter_gauss2D(movie3D, gausskernel)
%%  filter_gauss2D(movie3D, gausskernel) filters each frame independently with imgaussfilter
% INPUT: movie
% OUTPUT: filtered 3D movie

lastframe = size(movie3D, 3); 

    for fr = 1:lastframe
        output(:,:,fr) = imgaussfilt(movie3D(:,:,fr),gausskernel);
    end
end


%% Created: 13/02/21 (from Gent) 
% Updated: 