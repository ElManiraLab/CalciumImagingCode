% Confocal Video

% Load VSDmov.data. Use just one movie
Confocal = VSDmov.data(:,:,:,1);

v = VideoWriter('confocal.avi');
v.FrameRate = 9.5238;
open(v);

for i = 1:200
    imagesc(Confocal(:,:,i)); colormap('bone')
frame = getframe;
   writeVideo(v,frame);
end 

close(v);

% Load tail movement video

v = VideoWriter('tail_movement.avi');
v.FrameRate = 37.5714;
v.Duration
open(v);


for i = 52:840 % all videos should have shame length
    imagesc(frames_laser(:,:,i)); colormap('bone')
frame = getframe;
   writeVideo(v,frame);
end 

close(v);
