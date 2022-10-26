function [outputArg1,outputArg2] = MaximumInRoi(VSDI,VSDmov)
%MAXIMUMINROI Summary of this function goes here
%   Detailed explanation goes here

a=['r','b','g','c','k','']
VSDI.crop.preview

IMG=VSDI.crop.mask.*VSDI.crop.preview;
imagesc(IMG); colormap('bone');
hold on
[row, col]=size(VSDI.condition)
for c=1
    idx=[]
idx=find(VSDI.nonanidx(:,c));
for idxcond=1:length(idx)
    A=[]
    A=VSDmov.data(:,:,:);%idx(idxcond));
    [zz xx yy] =size(A)
    B=[]
    for dd=1:yy
        B(:,:,dd)= A(:,:,dd).*VSDI.crop.mask;
    end
    Maxi = max(max(max(B)))
    [r,c,v] = ind2sub(size(B),find(B >Maxi*0.80));
%     pgon = polyshape(c,r)
%     plot(pgon)
    
p=plot(c,r, strcat(a( idxcond),'.'))
p.MarkerSize = 15
    hold on
end 
end 

end

% MAS FACIL 

SinBack=VSDmov.data(:,:,1:end-1);
SinBack=SinBack(:,:,[148]);
Maxi = max(max(SinBack))
[r,c,v] = ind2sub(size(SinBack),find(SinBack >Maxi*0.4));
p=plot(c,r, 'r.')


IMG=VSDI.crop.mask.*VSDI.crop.preview;
imagesc(IMG); colormap('bone');
hold on
SinBack=filtmov(:,:,1:end-1);
SinBack=SinBack(:,:, 68);
Maxi = max(max(SinBack))
[r,c,v] = ind2sub(size(SinBack),find(SinBack >Maxi*0.8));
p=plot(c,r, 'r.')
p.MarkerSize = 15
