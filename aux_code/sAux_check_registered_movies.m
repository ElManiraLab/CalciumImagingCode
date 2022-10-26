% TO CHECK WHETHER s01_importmov_basic_preprocess.m - point 01

for tri=1:VSDI.nonanidx
   imagesc(rawmov(:,:,1,tri)); colormap('bone');
   title(num2str(tri))
   pause
end

for tri=1:VSDI.nonanidx
   imagesc(registermov(:,:,1,tri)); colormap('bone');
   title(num2str(tri))
   pause
end