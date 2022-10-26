%% David Madrid . Last Rev 28/12/2022 

%% FOLDER ORGANIZATION

% The organization of the folder must be as followed. 
%------------------------------------------------------------------------------------%

%- GCaMP_Date

%-----Calcium_ourToolbox

%----------------------------Calcium_data ([path.data])
%-------------------------------------CZI --> save confocal/2P files form ZEN software (.czi)
%-------------------------------------dataCALCIUM
%-------------------------------------datamovies
%-------------------------------------datawaves ([path.waves])
%-------------------------------------Deeplabcut--> csv files of the coordinates 
%-------------------------------------TIFF --> Storage of the images save as image sequence (.tif) 
%-------------------------------------VIDEO_camera --> Video (.mp4)
%-------------------------------------VIDEO_deeplabcut --> Video after DEEPLABCUT labelling (.mp4)

%----------------------------rootpath
%------------------------------------- Code
%----------------------------------------------user_settings--> this file
%----------------------------------------------CALCIUM_Analysis--> Computational code
%----------------------------------------------Post_Analysis
%----------------------------------------------Other Scripts

%-----General_details_ Experiment--> (Optional).
%-----READ_ME--> (Optional).

%------------------------------------------------------------------------------------%

%% USER_SETTINGS pathways. 

path.Calcium_ourToolbox = pwd ;  % run the code from the inside of the rootpath folder

% Next they are define the folders where all the data will be storage.
% Change this if you change the directory of the GCaMP_Date folder

% path.data is the directory to the folder "Calcium_data"
path.data  = fullfile(path.Calcium_ourToolbox(1:end-8),'Calcium_data');
% path.waves is the directory to the folder "datawaves"
path.waves = fullfile(path.data,'\','datawaves');

%%
% ....................................................................................................................
% ....................................................................................................................

path.grouplist = path.data ;
addpath(genpath(path.Calcium_ourToolbox)) ;


%% General directories for the experimental data. If some of them are
% not in use you can comment them

% Tiff_directory = fullfile(path.data,'TIFF'); % directory for the TIFF folder
CSV_directory =  fullfile(path.data,'Deeplabcut'); % directory for the csv folder
% Videocsv_directory =  fullfile(path.data,'VIDEO_deeplabcut'); % directory for the labelled videos
Videoreal_directory = fullfile(path.data,'VIDEO_camera');  % directory for the videos
CZI_directory = fullfile(path.data,'CZI'); % directory for the czi folder


%% TESTING NUMBER OF FILES in the folder.
%If the number of files  match it will throught an error. 

a2=dir([CSV_directory '/*.csv']); csv_num=length(a2); 
% a3=dir([Videocsv_directory '/*.mp4']); videocsv_num=length(a3);
a4=dir([Videoreal_directory  '/*.mp4']); videoreal_num=length(a4);
a5=dir([CZI_directory '/*.czi']); czi_num=length(a5);

% equality_matrix = [czi_num csv_num videocsv_num videoreal_num ];
equality_matrix = [czi_num  csv_num videoreal_num ];

if all(equality_matrix == equality_matrix(1))
    disp('All folders have the same number of elements')
else
    error(['Error: Some folders have different number of elements. Check the number'...
        ' of elements in the folders CZI,' ...
        ' videocsv, videoreal and csv to find the problem']) 
end  

%% EXPERIMENT LIST CREATION. For this it will be use the name of the .czi files

A_cell = struct2cell(a5); % Spliting czi
A_cell = A_cell(1,:)';
list = strings(length(A_cell),1);

for i = 1:length(A_cell)
 x = split(A_cell(i,1),".czi");
 list{i,1} = x{1};
end 


% Condition columns. If you dont want to analyze some specific trial (NaN) or one
% of them has a specific condition (N = {2,3,...,n}). N = 1 is reserved for
% the standar condition. 
C1((1:length(list)),1) = 1; %Must be in order % @SET


%% TESTING THE NAME AND ORDER OF THE FILES in the csv folder,czi folder and video folder

A_cell = struct2cell(a2); % Spliting csv
A_cell = A_cell(1,:)';
list2 = strings(length(A_cell),1);

for i = 1:length(A_cell)
 x = split(A_cell(i,1),"DLC_");
 list2{i,1} = x{1}; % Creating list2 with the name of the csv
end 

A_cell = struct2cell(a4); %spliting video
A_cell = A_cell(1,:)';
list3 = strings(length(A_cell),1);

for i = 1:length(A_cell)
 x = split(A_cell(i,1),".mp4");
 list3{i,1} = x{1}; % Creating list3 with the name of the videos
end 


if all(list == list2) %&& all(list == list3)
    disp('All folders have the elements in order and with the same name')  
else
     error('Error: Some folders have elements in different order or with different names') 
end 

%% EXTARCTING NAMES

% % CZI
% 
% A_cell = struct2cell(a5); 
% A_cell = A_cell(1,:)';
% list1 = strings(length(A_cell),1);
% 
% for i = 1:length(A_cell)
%  list1(i,1) = A_cell(i,1); 
% end 

% CSV

A_cell = struct2cell(a2); 
A_cell = A_cell(1,:)';
list2 = strings(length(A_cell),1);

for i = 1:length(A_cell)
 list2(i,1) = A_cell(i,1); 
end 

% VIDEOS

A_cell = struct2cell(a4); 
A_cell = A_cell(1,:)';
list3 = strings(length(A_cell),1);

for i = 1:length(A_cell)
 list3(i,1) = A_cell(i,1); 
end 

% SAMPLING RATE






