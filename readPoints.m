[x,y] = getpts function pts = readPoints(image,Movie,Trial,n)
%readPoints   Read manually-defined points from image
%   POINTS = READPOINTS(IMAGE) displays the image in the current figure,
%   then records the position of each click of button 1 of the mouse in the
%   figure, and stops when another button is clicked. The track of points
%   is drawn as it goes along. The result is a 2 x NPOINTS matrix; each
%   column is [X; Y] for one point.
% 
%   POINTS = READPOINTS(IMAGE, N) reads up to N points only.

 I=(Movie.data(:,:,1,Trial));
if nargin <=3
    n = Inf;
    pts = zeros(2, 0);
else
    pts = zeros(2, n);
end
subplot(2,2,[1 3])
imagesc(I);     % display image
xold = 0;
yold = 0;
k = 0;
hold on;           % and keep it there while we plot
while 1
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    k = k + 1;
    pts(1,k) = xi;
    pts(2,k) = yi;
      if xold
          clear A
          plot(xi, [ yi], 'go-');  % draw as we go
                     xi=fix(xi)
      yi=fix(yi)
     
       A=Movie.data(yi,xi,:,Trial);
 A=A(1:end-1);
 A=squeeze(A);
 subplot(222)
 plot(A)
      else
          plot(xi, yi, 'go');         % first point on its own
           xi=fix(xi)
      yi=fix(yi)
     
       A=Movie.data(yi,xi,:,1);
 A=A(1:end-1);
 A=squeeze(A);
 subplot(222)
 plot(A)
          
      end
      if isequal(k, n)
          break
      end
      xold = xi;
      yold = yi;
      xi=fix(xi)
      yi=fix(yi)
      
 A=Movie.data(xi,yi,:,1);
 A=A(1:end-1);
 A=squeeze(A);
 plot(A)
  end