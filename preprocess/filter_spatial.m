function [output] = filter_spatial(input, P)
% SPATIAL AVERAGE FILTER OF THE IMAGE (2D input) OR OF THE FRAMES OF A
% MOVIE (3Dinput)

% performs "P x P" pixels average filter of 2D filter processing. We normally set the parameter "P" to 9.

% If The value defined as D(t, x, y) on the coordinates of (x, y) at the frame number "t" will be replaced with the averaged value - in the area "P x P".
% If input is an image (dim=2), returns an image filtered.
% If the input is a movie (3D), ut returns a movie in which each frame has
% been filtered independently

kernel = ones(P);
data_dim = length(size(squeeze(input)));

if data_dim ==2 
kernel = kernel / sum(kernel(:));
output = conv2(double(input), kernel, 'same');
end

if data_dim ==3
data3D = input; 
n_frame = size(input,3); 

for fr = 1:n_frame
frame = squeeze(data3D(:,:,fr));
kernel = kernel / sum(kernel(:));
localAverageImage = conv2(double(frame), kernel, 'same');
    
output(:,:,fr) = localAverageImage; 
end

end

end

function [output] = filter_Tcnst(input, weight)
% TIME-CONSTANT FILTER (as defined in brainvision) FOR TEMPORAL SMOOTHING
% OF WAVES (1Dinput) OR PIXELS OF MOVIES (3Dinput)
%
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
    x = size(input,1); y =size(input,2);
    npix= x*y;
    t= size(input,3);
    
    input_flat = reshape(input, npix, t);
    
    for px = 1:npix
    output_flat(px,:)= Tcnst1(input_flat(px,:), weight);
    end
    
    output = reshape(output_flat,x,y,t);
end

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
        outputT1(1) = (1-mult)*input(1)%+mult*outputT1(2);
        
    % time filter taking into account previous values
    for ii =  2:L
        outputT1(ii) = (1-mult)*input(ii)+mult*outputT1(ii-1);
    end

end