function [data3D,timeb]= BVimport(filepath, nframes)

% TRANSFORM csv FROM BRAINVISION ANALYZER INTO 3D DATA MATRIX

% filepath = 'C:\Users\User\Documents\UGent\data_extraction\code\fish_intensity\fish1h1int.csv';
%'nframes' - number of frames of activity (without the background frame)

% ROi HAVE TO BE extracted using saveEach from ROI  (and the csv 'wholewindowROI' stored in this
% folder)
% the difference in values between the extracted ROI and the image matrix
% is due to a systematic failure of BV when applying the ROI mask

% INPUT:
% 'filepath' of raw BV data
% number of frames in the movie
% OUTPUT:
    % data3D : 3D matrix {(x,y) image, (z)timeserie} in data3D(:,:1:680)
    % > background image is in data3D(:,:,682)
    % 'timeb' timebase 

lastframe = 2 + nframes; 

table = readmatrix(filepath, 'HeaderLines',7);
data3D = table(1:lastframe, 2:5134); %5134 is the number of pixels (that in the BV file are in a row)
Idata(1:nframes,:) = data3D(3:lastframe,1:end);

%the first row correspond to background value, 
% the first row (from row=2) to the Brain Vision timebase
backgr = table(1, 2:5134);
backgr= reshape(backgr(:),59,[]);
timeb = table (3:lastframe,1);
data3D = zeros(59,87,nframes);

for framei = 1:nframes
tempdata = squeeze(Idata(framei,:));
data3D(:,:,framei) = reshape(tempdata(:),59,[]);
end
data3D(:,:,lastframe) = backgr; 
end