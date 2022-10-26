function [movie_thres, alphachan] = movie_absthresh(movie,method,value)
%% TEST!!!! 

% INPUT: 
% 'movie' - 3D movie (with background in the last frame, that
% will be ignored)
% 'method'
% 'value' according to the method.
    % method = 'absvalue' >> value = [abs value] and all values with an
    % absolute value below will be zero-out
    % If 'method' = 'absmaxmin', value = [min max] 

% OUTPUT: 
% 'movie_thres' - zero-out values under threshold
% 'alphachan' - transparency logic matrix; 1~100% transparency
movie_thres = movie(:,:,1:end-1); %removes background

switch method
    case 'absvalue'

    movie_thres=movie; % preallocate
    
    movie_thres(abs(movie_thres) < value) =0;
    alphachan = movie_thres==0 ;
    
    case 'absmaxmin'
        
    case 'abspercent'     %here value is the percent (over 1) from the absolute value
    % Discard values above SD  !!!! 
    maxval = max(abs(movie(:)));
    movie_thres=movie; % preallocate
    
    movie_thres(abs(movie_thres) < maxval * value) =0;
    alphachan = movie_thres==0 ;
    
end
        movie_thres(:,:,end+1)= movie(:,:,end); % background back
end