location = 'C:\Users\Usuario\OneDrive\Escritorio\Calcium_ourToolbox\rootpath\TIFF\Experiment_18_09_21\1\*.tif';   
ds = imageDatastore(location)  

suma = 0

while hasdata(ds) 
    
    img = read(ds) ;             % read image from datastore
    % creates a new window for each image
   for i = 1
       suma = suma +1
   end 
   Hyperstack(:,:,suma) = img;
end

imtool(Hyperstack(:,:,80));

A = squeeze(Hyperstack(390,450,:));

plot(1:94,  A,'r')