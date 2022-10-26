light.likehood = round(CSV_file(:,10));
laserframeON = find( light.likehood,1,'first');
laserframeOFF = find( light.likehood,1,'last');

laserframe = laserframeON:laserframeOFF;
lasertime = time(laserframeON:laserframeOFF);

% Test 
% Read Video Object
vidObj = VideoReader(CALCIUM.Video_directory_loop);
% Create 4D a matrix with frames
framesON = read(vidObj,[1 laserframeON+1]);

% Use just one channel and go from 4D to 3D
frames_laserON = squeeze(framesON(:,:,2,:));
% Checking framer before, during and after laser
 subplot(231)
imagesc(frames_laserON(:,:,laserframeON-2))
title('Bfr laserON')

 subplot(232)
imagesc(frames_laserON(:,:,laserframeON))
title('laserframeON')

 subplot(233)
imagesc(frames_laserON(:,:,laserframeON+1))
title('Aft laserON')


time_sync = time(1:length(lasertime));

CALCIUMroiTS.deeplabcut.light.likehood= light.likehood;
CALCIUMroiTS.deeplabcut.laserframe = laserframe; 
CALCIUMroiTS.deeplabcut.lasertime = lasertime;
CALCIUMroiTS.deeplabcut.laserframeON  = laserframeON ;
CALCIUMroiTS.deeplabcut.laserframeOFF  = laserframeOFF ;
CALCIUMroiTS.deeplabcut.lasertimesync = time_sync; %
CALCIUMroiTS.deeplabcut.tail1.x = CSV_file(laserframe,2);
CALCIUMroiTS.deeplabcut.tail1.y = CSV_file(laserframe,3);
CALCIUMroiTS.deeplabcut.tail2.x = CSV_file(laserframe,5);
CALCIUMroiTS.deeplabcut.tail2.y = CSV_file(laserframe,6);

CALCIUMimg('savewave', CALCIUMroiTS, CALCIUMroiTS.ref,list,nfish);

figure(2)
difff = CALCIUMroiTS.deeplabcut.tail1.x (1) - CALCIUMroiTS.deeplabcut.tail2.x(1);
plot(time_sync,CALCIUMroiTS.deeplabcut.tail1.x,'b'), hold on
plot(time_sync,CALCIUMroiTS.deeplabcut.tail2.x +difff,'r') 
title(char(CALCIUM.list(nfish,1)))
ylabel('\Deltax Coordinates tail (red) + arduino(blue) [px]')
xlabel ('Time [s]')
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Tail_Movement_',char(CALCIUM.list(nfish,1)))))
close

x = CALCIUMroiTS.deeplabcut.lasertimesync;
y1 = CALCIUMroiTS.deeplabcut.tail1.x;
y2 = CALCIUMroiTS.deeplabcut.tail2.x +difff;


for i = 1:length(CALCIUMroiTS.roi.labels)
yyaxis left
plot(x,y1,'g','LineWidth',3); hold on
plot(x,y2,'r','LineWidth',3); hold off

r = -waves(:,i,1);
timev = CALCIUM.timebase;
yyaxis right
plot(timev,r,'b','LineWidth',3);

yyaxis left
title(strcat('Calcium activity-Tail movement','-',CALCIUM.ref));
xlabel('Time [s]')
ylabel('\Deltax Coordinates tail3 [px]')

yyaxis right
ylabel('% \Delta F')
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),...
    strcat(char(CALCIUMroiTS.roi.labels(i)),'_',char(CALCIUM.list(nfish,1)))))

legend('arduino swimming','spont swimming', 'GCaMP6f')

end 
close

CALCIUM.delay = x(end)-timev(end);
CALCIUMimg('save',CALCIUM,[],list,nfish);

clear 

