%% RASTER. 
%
% For a specific colormap.
% For finding a specific color go to https://htmlcolorcodes.com/es/

% For lateral cells
lengthH = 1000;
neutral = [243, 220, 218]./255;
blue = [87, 167, 179]./255;
colors_p = [linspace(neutral(1),blue(1),lengthH)', linspace(neutral(2),blue(2),lengthH)',...
    linspace(neutral(3),blue(3),lengthH)'];
% For medial cells
lengthH = 1000;
neutral = [243, 220, 218]./255;
violet = [174, 74, 132]./255;
colors_p2 = [linspace(neutral(1),violet(1),lengthH)', linspace(neutral(2),violet(2),lengthH)',...
    linspace(neutral(3),violet(3),lengthH)'];


h = imagesc(perloc); colormap(colors_p); colorbar ('southoutside')
 
h2 = imagesc(perloc); colormap(colors_p2);colorbar ('southoutside')
 

%% CURVE WITH STANDAR DEVIATION SHADOW
%medial 

k = -CALCIUMroiTS.diff_perc03.data(:,4);
j = -CALCIUMroiTS.diff_perc03.data(:,5);
P= -CALCIUMroiTS.diff_perc03.data(:,6);
l= -CALCIUMroiTS.diff_perc03.data(:,7);
e= -CALCIUMroiTS.diff_perc03.data(:,5);
g= -CALCIUMroiTS.diff_perc03.data(:,6); 
for i=1:length(k)
%avg1(i)=(k(i)+j(i))/2;
% avg1(i) = (k(i)+j(i)+P(i)+l(i)+e(i)+g(i))/6;
avg1(i) = (k(i)+j(i)+P(i)+l(i))/4;
end


y = avg1; % your mean vector;
x = CALCIUMroiTS.diff_perc03.times;
for i=1:length(avg1)
    stdv(i) = std([k(i),j(i)])
end 
curve1 = y + stdv;
curve2 = y - stdv;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [236, 232, 235]./255);
hold on;
plot(x, y,'color', violet, 'LineWidth', 2);


% laterals
 figure

% o = -CALCIUMroiTS.diff_perc03.data(:,7);
s = -CALCIUMroiTS.diff_perc03.data(:,8);
d = -CALCIUMroiTS.diff_perc03.data(:,9);
k = -CALCIUMroiTS.diff_perc03.data(:,10);
j = -CALCIUMroiTS.diff_perc03.data(:,11);
P= -CALCIUMroiTS.diff_perc03.data(:,12);
l= -CALCIUMroiTS.diff_perc03.data(:,13);


for i=1:length(k)
% avg2(i)=(o(i)+s(i)+d(i))/3;
avg2(i) = (k(i)+j(i)+P(i)+l(i)+s(i)+d(i))/6;
end

y = avg2; % your mean vector;
x = CALCIUMroiTS.diff_perc03.times;
for i=1:length(avg2)
    stdv(i) = std([o(i),s(i),d(i)]);
end 
curve1 = y + stdv;
curve2 = y - stdv;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [236, 232, 235]./255);
hold on;
plot(x, y,'color', blue, 'LineWidth', 2);


%% Amplitude plot

% take values of amplitude for one episode of swimming and transformed in a
% vector 
episode1 = [];
episode2= [];
episode3= [];
episode4= [];


boxplot(episode1)
ylim([10 70])
boxplot(episode2)
ylim([10 70])

boxplot(episode3)
ylim([10 70])

boxplot(episode4)
ylim([10 70])
%% BOX PLOT


