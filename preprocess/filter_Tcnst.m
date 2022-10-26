function [output] = filter_Tcnst(input, weight)
% [output] = filter_Tcnst(input, weight)
% TIME-CONSTANT FILTER (as defined in brainvision) FOR TEMPORAL SMOOTHING
% OF WAVES (1Dinput) OR PIXELS OF MOVIES (3Dinput)
% INPUT: wave or movie to filter (if movie, the background is in the last
% frame and will be excluded from the computation)

% The function that does the actual calculation is the attached  Tcnst1. 
% It uses 'weight' number of previous timepoints to smooth. The first points are smoothed in a reversed way, using next timepoints 
% input: one timeserie( Wave)
% weight: of the Tcnst filtering ~ timespan used to smooth (see below)
% output: filtered timeserie or movie

% BRAINVISION DEFINITION:
% Change the Temporal Information of the Wave Chart
% 
% (D(t): value of wave chart before calculation
%  D'(t): value of wave chart after calculationt: time
% 
% When Tcnst=1, D'(t) = D(t)*1
% When Tcnst=2, D'(t) = D(t)*0.9 + D’(t-1)*0.1
% When Tcnst=3, D'(t) = D(t)*0.8 + D’(t-1)*0.2
% When Tcnst=4, D'(t) = D(t)*0.7 + D’(t-1)*0.3
% ???
% When Tcnst=10 D'(t) = D(t)*0.1 + D’(t-1)*0.9

if length(size(input))== 2 % waves, although being 1D, they are represented in matlab as 2-D (to indicate if it's a row or a column)
    [output] = Tcnst1(input, weight);
end

if length(size(input))== 3 % If it is a movie, the filter is applied to each pixel
    movie = input(:,:,1:end-1); %leave the backgrond out of the computation
    
    x = size(movie,1); y =size(movie,2);
    npix= x*y;
    t= size(movie,3);
    
    input_flat = reshape(movie, npix, t);
    
    for px = 1:npix
    output_flat(px,:)= Tcnst1(input_flat(px,:), weight);
    end
    
    output = reshape(output_flat,x,y,t);
    output(:,:,end+1) = input(:,:,end); % restore the background frame

end %end if movie

end

function [outputT1] = Tcnst1(input, weight)
    
    if weight ==1
    mult = 1;
    end

    if weight > 1 
       mult = (weight - 1) /10;
    end

    %When Tcnst=10 D'(t) = D(t)*0.1 + D’(t-1)*0.9
    L = length(input);
    outputT1 = NaN(L,1); 
    
    % fix first point, a reverse filter is applied
        outputT1(1) = (1-mult)*input(1);%+mult*outputT1(2);
        
    % time filter taking into account previous values
    for ii =  2:L
        outputT1(ii) = (1-mult)*input(ii)+mult*outputT1(ii-1);
    end

end

%% Created: 13/02/21 (from old code)
% Updated: 17/02/21 to leave the background out the computations