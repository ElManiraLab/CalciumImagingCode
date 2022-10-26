%% UPLOADING Z-STACK
user_settings

data = bfopen(char("C:\Users\Usuario\OneDrive\Escritorio\TEMP_DATA\GCaMP_220419\Calcium_ourToolbox\Calcium_data\ZSTACK\Experiment-201.czi"));

for j = 1:length(data{1,1})
  zstack (:,:,j) = data{1,1}{j,1};
end
 zstack2 = im2double(zstack);

    nfish =4; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;




  %% ERROR
 for i=1:size(zstack2 ,3)
 err(i) = immse(zstack2(:,:,i),CALCIUM.mean);
  end 
  find( err==min( err))


%% MULTIPLICATION
for i=1:54
 multi(:,:,i) = CALCIUM.mean.*zstack2(:,:,i);
end 

index = squeeze(sum(sum(multi)));
find(index==max(index))

%% SUBSTRACTION
for i=1:54
 substrac(:,:,i) = abs(CALCIUM.mean-zstack2(:,:,i));
end 

index = squeeze(sum(sum(substrac)));
find(index==min(index))


%%
 subplot(121)
 imagesc( zstack2(:,:,34))
 subplot(122)
 imagesc(CALCIUMmov.data(:,:,1))

 for i=1:size(zstack2 ,3)
      imagesc( zstack2(:,:,i))
      pause(0.2)
      
 end 


for i=1:size(zstack2 ,3)
 ssimval(i) = ssim(zstack2(:,:,i),CALCIUM.mean);
end 

find(ssimval==max(ssimval))




