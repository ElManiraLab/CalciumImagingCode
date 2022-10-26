%% AUTOMATIC LASER ON/OFF DETECTION
% Read Video Object
vidObj = VideoReader('C:\Users\Usuario\OneDrive\Escritorio\video\211029_CB_2DLC_resnet50_GCaMPOct29shuffle1_1030000_filtered_labeled.mp4');
% Create 4D a matrix with frames
frames = read(vidObj,[1 vidObj.NumFrames]);
% Use just one channel. It seems like with channel 2 the  difference are bigger
frames_bw = frames(:,:,2,:);
% Go from 4D to 3D
frames_laser = squeeze(frames_bw);
% Average first seconds without light
nolaser = mean(mean(frames_laser(:,:,1))); %@ SET FRAMES (ei. 100)
% Loop to find the first frame with the laser comparing
% pixel intensity with the no laser average
for i = 1:vidObj.NumFrames
   if  mean(mean(frames_laser(:,:,i))) < 5 %@ SET THRESHOLD
       laserframe(i)= 0;
   elseif  mean(mean(frames_laser(:,:,i))) > 5%@ SET THRESHOLD
       laserframe(i)= 1;
   end 
end 

% Finding first non zero frame
 laserframeON = find( laserframe,1,'first');
 laserframeOFF = find( laserframe,1,'last');
 
 % Checking framer before, dufring and after laser
 subplot(141)
imagesc(frames_laser(:,:, laserframeON-1))
title('Bfr laser')

 subplot(142)
imagesc(frames_laser(:,:, laserframeON))
title('laserframeON')

 subplot(143)
imagesc(frames_laser(:,:, laserframeON +1))
title('Aft laser')

 subplot(144)
imagesc(frames_laser(:,:, laserframeOFF))
title('No Laser')

display(laserframeON)

%% SAMPLING RATE CONVERSION
% Sampling rate of the video (phone/camera)
fps = 60; 
% Interval between frames (ms)
ibf = 1000/ fps ;
% Timebase vector
x = 0:16.667:(vidObj.NumFrames*16.667)-1;
 



