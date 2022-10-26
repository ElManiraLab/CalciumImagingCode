%% David Madrid . Last Rev 14/04/2022

%%                                 CLEAR AND CLOSE ALL

clear
clc
close all
user_settings;
list


temp_table(:,1) = table(list); % First column of the table --> list
temp_vector = 1:length(list); % Second column of the table --> 1:length(list)
temp_table(:,2) = table(temp_vector'); 
temp_table.Properties.VariableNames = {'List','nfish'};  
disp(temp_table) % Display the table as info

%%                              LOADING FISH

for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           SORTING NEURONS LEFT TO RIGTH, MID AND LATERAL REFERENCE AND NORMALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Mean x cordinate and y cordinate (center) per each cell. For each ROI
    % created in the CALCIUM_Analysis script, we will calculate the center in
    % the x-y plane (coordinates). the variable manual_poly contains the
    % coordinates for each point the roi have, so averaging the x valus and the
    % y values will give us the center of the cell

    xpos_N = NaN(1,size(CALCIUM.roi.manual_poly,1)); % Prelocate the vector for the mean x values
    ypos_N = NaN(1,size(CALCIUM.roi.manual_poly,1)); % Prelocate the vector for the mean y values

    for i = 1:size(CALCIUM.roi.manual_poly,1) % the size will be equal to the number of neurons
        temp = CALCIUM.roi.manual_poly{i}; % for each x and y coordinate of one neuron
        xpos_N(i) = mean (temp(:,1)); % the first column is the x coordinates of each point
        % of one ROI of one neuron
        ypos_N(i)= mean (temp(:,2)); % the second column is the y coordinates of each point
        % of one ROI of one neuron
        clear temp
    end
    num_neurons = xpos_N; % vector containing the average x position for each neuron
    y_neurons = ypos_N; % vector containing the average y position for each neuron

    sorted = sort(num_neurons); % sort the x position (pixel) of
    % the neurons (ascendent order: smallest to highest): Ex-> (30px 2px 10px) = (2px 10px 30px)

    clear xpos_N


    if isfield(CALCIUM,'Midreference') % checking the condition of the existence of a mid reference

        % Because the midline can be a lil bit tilted, we can not take the
        % average value of the x colum of the midline. The x value that we take
        % must be releated with the closest y posistion in the midline to the
        % cell. i.e. if the the midline is tilted some degrees to the leftnand
        % just take the x averages we can have cells on the top part of the
        % picture that are closer to that point and cells n the bottom that are
        % futher away, when in realiy all of then are at same distance of their
        % respective point of the midline. Because the x poisition of the
        % midline is varying with the y cordinate we have to select a cell,
        % find the y cordinate that is closer in the midline, and them used the
        % x value of that coordinate like the midline for that cell. If the
        % midline is completely vertical, doenst matter all the y cordinates of
        % the midline have the same x coordinate.


        for i = 1:length(CALCIUM.Midreference.manual_poly) % carefull now
            % is made just for one ROI (MID REFERENCE), the midline
            temp = CALCIUM.Midreference.manual_poly{i}; % will contain the x and y values of all the
            %point that contains the midline ROI. This temp variable will be
            %use later. dont erase it.
        end

        %Prelocate
        midline=NaN(1,length(ypos_N)); % Will be the x value that we will use
        % as a midline per each cell

        for i = 1:length(ypos_N) % Loop across cells
            temp1 = abs(temp(:,2) - ypos_N(i)); % we find the y coordinate in the midline that is closer to cell
            index = temp1==min(temp1);
            index = find(index);
            midline(i) = temp(index(1),1); % x value, that is why we use the fist colum of temp
        end
        clear temp1 index


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT TO RIGHT ORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        index_order = NaN(1,length(sorted));

        for i=1:length(sorted)
            index_order(i)= find(num_neurons==sorted (i)); % find the order
            %in the actual variable: Ex-> find the index for the smallest value in
            %the variable and tell me which index is it, plus save that index like
            %the first one in the new variable index_order. So the fist value of
            %this new variable will be the index of the sammlest value of the
            %original value. If the fist value of index_order is 3, means that the
            %smallest value in this variable num_neurons is in the
            %position number 3 when is not sorted
        end


        perloc = zeros(size(CALCIUMroiTS.diff_perc03.data,1),...
            size(CALCIUMroiTS.diff_perc03.data,2)); %preloc the new matrix
        name_perloc = cell(1,length(CALCIUM.roi.labels)); % Preloc for organazing
        % also the name of the neurons

        for i=1:length(sorted)
            perloc(:,i) = CALCIUMroiTS.diff_perc03.data(:,index_order(i)); % the row column
            % will be the row with the smallest  x position and the last row will be
            % the rox with the highest x positon
            name_perloc(i) = CALCIUM.roi.labels(index_order(i));
        end


        % perloc = flipud(perloc');
        perloc=-perloc'; % The column is the one with the lower x position of the pixel
        name_perloc = name_perloc';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBSTRACTING MID REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % For normalizing to the midline, we substract first the correspondient x
        % value of the midline to each cell

        distance_to_mid = NaN(1,length(name_perloc));

        for i=1:length(name_perloc)

            distance_to_mid (i) = num_neurons(index_order(i)) - midline(index_order(i));% This
            % is the distance to the midline already sorted from left to rigth; This is the distance (x distance)
            % to the midline for each cell, sorted from left to right. So
            % negative values means that are in the left hemisphere, and
            % positive values means that are in the right hemisphere. Bigger
            % the absolute value bigger the distance to the midline. the value
            % 0 corresponds to the midline. this is normalize yet to the most
            % lateral cell. As you can see the number are sorted form smallest
            % (most negative) to biggest (most positive). that means the cells
            % are sorted from left to right according to their possition in the
            % x line. THIS VARIABLES HAVE BEEN APLIED ALREADY THE INDEX ORDER
            % SO IS ALREADY ORDER.

        end

    end


    %%%%%%%%%%%%%%% NORMALIZING LATERAL REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if   isfield(CALCIUM,'Latreference') % Checking for existence of a lateral reference

        for i=1:length(CALCIUM.Latreference.manual_poly) % carefull now
            % is made just for un ROI (LATERAL REFERENCE)
            temp2 = CALCIUM.Latreference.manual_poly{i}; % first column: x coordinates of the
            %lateral reference, second column: y coordinates of the lateral
            %reference
        end

        Lateral_position_X = mean(temp2(:,1)); % this is the average x position of an
        % identificable lateral reference.
        Lateral_position_Y = mean(temp2(:,2)); % this is the average y position of an
        % identificable lateral reference.


        % now we substract this distance to the midline. we will find first
        % which y coordinate of the midline is closer to the y coordinate of
        % the lateral reference, the same thing we did to find the midline for
        % all the cells above

        temp1 = abs(temp(:,2) - Lateral_position_Y); % we find the y coordinate in the midline that is closer
        % to the y coordinate of cell
        index = find(temp1==min(temp1)); % index of the closer match
        midline_for_lateral = temp(index,1); % x value for the midline for the lateral reference

        lateral_norm = abs(midline_for_lateral-Lateral_position_X); % the distance to the midline that will be interpret as 1.

    end


   

    

    % Last step is divide the distance_to_mid by lateral_norm   and  all the
    % distance to the midline will be normalized to a reference distance which
    % is the distance to the midline of an identificable neuron.
    %     -> 0 will be the midline
    %     -> 1 will be the distance to the identificable neuron

    CALCIUM.Normalized_x_distance = distance_to_mid./lateral_norm;% this value will no be in order if the midline is not completly vertical
    CALCIUM.index_order =  index_order;
    CALCIUM.perloc = perloc; 
    CALCIUM.name_perloc = name_perloc; 


     % To find the position of the midline regarding to the neurons 

     if isempty(find(CALCIUM.Normalized_x_distance== min(abs(CALCIUM.Normalized_x_distance))))
         temp = find(CALCIUM.Normalized_x_distance== -min(abs(CALCIUM.Normalized_x_distance))); 
         CALCIUM.midline_respect_neurons = (temp+temp+1)/2; 
     else
         find(CALCIUM.Normalized_x_distance== min(abs(CALCIUM.Normalized_x_distance)))
          CALCIUM.midline_respect_neurons = (temp+temp-1)/2; 
     end
   

    CALCIUMimg('save',CALCIUM,[],list,nfish);

    % This analyses allow us to have values of Normalized_x_distance greater
    % than 1 because we can have cells with distances to the mid lines farther
    % away than the identificable neuron.


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                                    CALCIUM IMAGING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %%
%     % FILTERING CALCIUM IMAGING (MEDIAN FILTERING)
%     %CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
%     
%     CALCIUMroiTS.diff_perc03.data_No_filtered = CALCIUMroiTS.diff_perc03.data; %This will
%     % contain the data that is not filtered
% 
%     figure('units','normalized','outerposition',[0 0 1 1])
%     subplot(211)
%     plot(-CALCIUMroiTS.diff_perc03.data)
%    
% 
% 
% 
%     prompt = {'Do you wanna filter the data?:'};
%     dlgtitle = 'Input';
%     dims = [1 35];
%     definput ={'Y'};
%     answer = inputdlg(prompt,dlgtitle,dims,definput);
% 
%     if isequal(answer,{'Y'})
%         median = 5;
%         for i = 1:size(CALCIUMroiTS.diff_perc03.data,2)
%             temporal_wave(:,i) = medfilt1(CALCIUMroiTS.diff_perc03.data(:,i),median);
%         end
%         subplot(212)
%         plot(-temporal_wave)
%         
%         prompt = 'If you want to keep filtering type "Y", if you are satisfy type "N" : ';
%         str = input(prompt,'s');
% 
%         if isempty(str)
%             str = 'N'; % you are done
%             close
%         else
%             while isequal(str,'Y')
% 
%                 prompt = {'Change median threshold?'};
%                 dlgtitle = 'Input';
%                 dims = [1 35];
% 
%                 definput ={sprintf('%.0f',median)};
%                 answer = inputdlg(prompt,dlgtitle,dims,definput);
%                 median = str2double(answer(1));
% 
%                 for i = 1:size(CALCIUMroiTS.diff_perc03.data,2)
%                     temporal_wave(:,i) = medfilt1(CALCIUMroiTS.diff_perc03.data(:,i),median);
%                 end
% 
%                 subplot(212)
%                 plot(-temporal_wave)
% 
% 
%                 prompt = 'If you want to keep filtering type "Y", if you are satisfy type "N" : ';
%                 str = input(prompt,'s');
%             end
% 
%         end
%         CALCIUMroiTS.diff_perc03.data = temporal_wave;
%         CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish); %Saving the data filtered
% 
%     else
%     end
%     close all
%     clear temporal_wave

    %%                         AMPLITUDE OF THE CALCIUM IMAGING

    CALCIUMpeak = 60; % @SET value

    for i = 1:size(CALCIUMroiTS.diff_perc03.data,2) % Loop across all the cells

        % First we indentify the peaks of the CALCIUM WAVE
        figure('units','normalized','outerposition',[0 0 1 1])
        CALCIUMpeak = 60; % @SET value
        findpeaks(-CALCIUMroiTS.diff_perc03.data(:,i),'MinPeakProminence',CALCIUMpeak); hold on

        %Maybe is posible that you want to change the value of the threshold
        %for the peak per each cell so we will open a prompt for introducing
        %the new threshold value

        prompt = {'Enter Threshold:','Are you satisfy with the result?:'};
        dlgtitle = 'Input';
        dims = [1 35];

        definput ={sprintf('%.0f',CALCIUMpeak),'N'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        CALCIUMpeak = str2double(answer(1));
        cla

        % Conditional to repit again the process if the answer
        % is no. In this part, you can change dynamicly the value of the
        % threshold in the GUI

        while isequal(answer(2),{'N'})
            findpeaks(-CALCIUMroiTS.diff_perc03.data(:,i),'MinPeakProminence',CALCIUMpeak); hold on
            prompt = {'Enter Threshold:','Are you satisfy with the result?:'};
            dlgtitle = 'Input';
            dims = [1 35];

            definput ={sprintf('%.0f',CALCIUMpeak),'N'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            CALCIUMpeak = str2double(answer(1));
            cla
        end

        findpeaks(-CALCIUMroiTS.diff_perc03.data(:,i),'MinPeakProminence',CALCIUMpeak); hold on




        % Once the threshold is set, we can remove some peaks. First we create
        % the variable peaks and locs
        [pksCAL,locsCAL] = findpeaks(-CALCIUMroiTS.diff_perc03.data(:,i),'MinPeakProminence',CALCIUMpeak);

        %Creating a dialog box



        prompt = {'You want to remove some peaks:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput ={'Y'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        while isequal(answer(1),{'Y'})
            [x_time,~] = getpts;

            % Finding closest match
            removing = NaN(1,length(x_time)); % Prelocate
            for j=1:length(x_time)
                removing(j) = find(abs(locsCAL-x_time(j))==min(abs(locsCAL-x_time(j))));
            end

            % Removing this index values from the peaks and locs vectors
            pksCAL(removing) = [];
            locsCAL(removing) = [];
            cla
            findpeaks(-CALCIUMroiTS.diff_perc03.data(:,i),'MinPeakProminence',CALCIUMpeak);
            scatter(locsCAL,pksCAL,30,'filled','g')
            pause(0.5)
            prompt = {'You want to remove some peaks:'};
            dlgtitle = 'Input';
            dims = [1 35];
            definput ={'Y'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
        end

        % Testing

        scatter(locsCAL,pksCAL,30,'filled','g')

        clear x_time y_time removing
        clear x_time

        % Because each cell can have diferent number of peaks, we will store
        % the data in a structure


        my_field = strcat('CELL',num2str(i)); % Creating dinamicly the

        AMPLITUDE_Y.(my_field) = pksCAL; % we put
        AMPLITUDE_X.(my_field) = CALCIUMroiTS.diff_perc03.times(locsCAL); % we put
        % Finally we clear the pksCAL and locsCAL

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%% Finding onset of the calcium wave
        for kk = 1:length(pksCAL) % Loop across waves (peaks)

            % FINDING ONSET REFERING TO THE A % OF THE PEAK
            %{
      rising_threshold = 0.2; % Amplitud must be 20% of the original peak
      rising_points = pksCAL(kk)*0.2; % value of the amplitude of the rising point of the first wave (kk) 
      for the first neuron (i)
            %}


            % Because there can be several waves, there will be several rising_points
            % that will be equal to the rising_points of one cell. To avoid this, we
            % have to restrict the time points where we will find this rising_points to
            % points that just precede the peak that we are analysing

            time_segment(1) = locsCAL(kk)*CALCIUM.srate-20*CALCIUM.srate; % start of the segment
            % ntimes*sampling rate before the peak
            time_segment(2) = locsCAL(kk)*CALCIUM.srate; % end of the segment is the
            % moment in seconds of the peak

            new_vector_time = time_segment(1):CALCIUM.srate:time_segment(2);
            for kkk = 1:length(new_vector_time)
                index_segment (kkk)= find(round(CALCIUMroiTS.diff_perc03.times,4)==...
                    round(new_vector_time(kkk),4));  % indexs of the time CALCIUMroiTS.diff_perc03.times
                % that correspond with ntimes*sampling rate before the peak. we will
                % use this to index the CALCIUMroiTS.diff_perc03.data
            end


            % FINDING ONSET REFERING TO THE A % OF THE PEAK
            %{
      temp_rising = abs(rising_points-abs(CALCIUMroiTS.diff_perc03.data(index_segment,i))); 
      index_rising = find(temp_rising==min(temp_rising)); 
      time_rising = new_vector_time(index_rising(end)); % these will be the temporal indexs
            %}

            segment_data = CALCIUMroiTS.diff_perc03.data(index_segment,i);

            for jj = 1:length(CALCIUMroiTS.diff_perc03.data(index_segment,i))-1% rate of change
                slope(jj)=segment_data(jj)...
                    -segment_data(jj+1);
            end

            thresholh_segment = 6;
            sumation = 1;
            maxim_threshold = 0;
            index = min(find(slope>thresholh_segment));
            while size(index,2) == 0 % Can be that the original thershold was too big and there is not
                % any nuber above it. It will create a problem in the next part
                % of the code. For this reason
                thresholh_segment = thresholh_segment-sumation;
                index = min(find(slope>thresholh_segment));
                maxim_threshold = 1;


                % This while loop wil continue until there is an aceptable value
                % for the threshold

            end
            if maxim_threshold==1 % To inform to the user
                msgbox('Highest possible value for threshold reached')
                pause(2)
                clear maxim_threshold

            end

            CALCIUM.ONSET_X.(my_field)(kk) = find(round(CALCIUMroiTS.diff_perc03.times,4)==...
                round(new_vector_time(index),4));

            scatter(CALCIUM.ONSET_X.(my_field)(kk),-segment_data(index),30,'filled','r')

            prompt = {'Do you wanna set a new threshold','Which one?','Do you prefer using the value of the peak?'};
            dlgtitle = 'Input';
            dims = [1 35];

            definput ={'N',num2str(thresholh_segment),'N'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            answer1 =  answer(1);
            answer3 = answer(3);

            if isequal(answer3,{'Y'}) %if your wave comes after a second one, maybe
                % you dont have a proper onset, in this case, take the peak
                CALCIUM.ONSET_X.(my_field)(kk) = locsCAL(kk);
                scatter(CALCIUM.ONSET_X.(my_field)(kk),pksCAL(kk),30,'filled','b')

            elseif isequal(answer1,{'Y'})
                while isequal(answer1,{'Y'})
                    thresholh_segment = str2double(cell2mat(answer(2)));
                    index = min(find(slope>thresholh_segment))
                    maxim_threshold = 0;
                    while size(index,2) == 0 % Can be that the original thershold was too big and there is not
                        % any nuber above it. It will create a problem in the next part
                        % of the code. For this reason
                        thresholh_segment = thresholh_segment-sumation;
                        index = min(find(slope>thresholh_segment));
                        maxim_threshold = 1;


                        % This while loop wil continue until there is an aceptable value
                        % for the threshold

                    end
                    if maxim_threshold==1 % To inform to the user
                        msgbox('Highest possible value for threshold reached')
                        pause(2)
                        clear maxim_threshold

                    end

                    CALCIUM.ONSET_X.(my_field)(kk) = (find(round(CALCIUMroiTS.diff_perc03.times,4)==...
                        round(new_vector_time(index),4)));

                    scatter(CALCIUM.ONSET_X.(my_field)(kk),-segment_data(index),30,'filled','b')

                    prompt = {'Do you wanna set a new threshold','Which one?'};
                    dlgtitle = 'Input';
                    dims = [1 35];

                    definput ={'Y',num2str(thresholh_segment)}; 
                    answer = inputdlg(prompt,dlgtitle,dims,definput);
                    answer1 =  answer(1);

                end
            end

            CALCIUM.ONSET_X.(my_field)(kk) = CALCIUM.ONSET_X.(my_field)(kk)*CALCIUM.srate;

            saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
                char(CALCIUM.list(nfish,1)),strcat('Peaks_and_onset',my_field, char(CALCIUM.list(nfish,1)))))

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clear locsCAL pksCAL sumation
        w = waitforbuttonpress; % Once the scatter is plotted just click any key to go to next celL
        close all


    end

    %Saving

    CALCIUM.AMPLITUDE_Y = AMPLITUDE_Y;
    CALCIUM.AMPLITUDE_X = AMPLITUDE_X;
    CALCIUMimg('save',CALCIUM,[],list,nfish);

    %%                   EXTRACTING THE BASLINE FOR THE CALCIUM IMAGING

    % Extracting the baseline for each cell in order to calculate the correct
    % amplitude to the peak.

    %Prelocation of the variables
    baseline_frames = NaN(2,size(CALCIUMroiTS.diff_perc03.data,2));
    baseline_CALCIUM = NaN(1,size(CALCIUMroiTS.diff_perc03.data,2));
    xi_baseline_cutted = NaN(2,size(CALCIUMroiTS.diff_perc03.data,2));

    % Loop across cells
    for i = 1:size(CALCIUMroiTS.diff_perc03.data,2)
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,i))
        disp('Click once at the begining of the baseline, and twice at the end')
        [x_time,~] = getpts; %@

        % Control just to be sure we select 2 points
        if length(x_time)>2
            error('You click more the two points')
        elseif length(x_time)<2
            error('You click just the first point')
        end


        baseline_frames (:,i) = x_time; %Saving value in the variable baseline_frames
        %
        %     round(baseline_frames) % to have the value of the frame

        for j=1:size(baseline_frames,1) %Finding the index in the time vector that is
            % closer to the values we just pick

            substraction = abs(CALCIUMroiTS.diff_perc03.times-baseline_frames (j,i));
            xi_baseline_cutted (j,i) = find(substraction==min(substraction));
            clear substraction % Erase temporal variable
        end

        %Computing value of the baseline

        %temp will be a temporary vector that include the values of the calcium
        %wave btw the specific index that we selected. we will calculate the
        %average of this to have the avrg baseline for EACH cell.

        temp = mean(-CALCIUMroiTS.diff_perc03.data([xi_baseline_cutted(1,i):xi_baseline_cutted(2,i)],i));
        baseline_CALCIUM (i) = temp;
        clear temp

        close
    end

    user_settings
    CALCIUM.baseline_amplitude = baseline_CALCIUM;
    CALCIUMimg('save',CALCIUM,[],list,nfish);  % If you want to save the raw movie


    %%                    BASELINE TO PEAK FOR EACH CELL CALCIUM IMAGING

    %Having the value of the baseline for each cell and having also all the
    %peaks in amplitude for each cell, we can calculate now the
    %baseline-to-peak distance for each cell. the number of distance per cell
    %will be equal to the number of amplitude peaks. Because each cell can have
    %different peaks we can not save all the cells in a matrix, we have to use
    %a structure

    % Loop across cells
    for k = 1:size(CALCIUMroiTS.diff_perc03.data,2)

        my_field = strcat('CELL',num2str(k)); % Creating dinamicly the name of the
        % cell

        % Performing the baseline to peak per each cell, for that we use the
        % variable that we created perviously called AMPLITUDE_Y.(my_field)
        % where my field is CELL1,CELL2,... So we are performing the
        % substraction for the baseline of each cell with all the peaks of the
        % same cell to have the baseline to peak for that cell, and we are
        % looping across cell

        BASELINE_TO_PEAK.(my_field)  = CALCIUM.AMPLITUDE_Y.(my_field)-CALCIUM.baseline_amplitude(k);

    end

    % Saving

    CALCIUM.BASELINE_TO_PEAK = BASELINE_TO_PEAK;
    CALCIUMimg('save',CALCIUM,[],list,nfish);

    %% SET UP A THRESHOLD FOR THE CALCIUM (IN THIS RECORDING)
    % When you are detecting the peaks of each cell the y axis of each cell
    % is diferent and then you dont know if the peak that you are selecting
    % is big or small. Set a threshold that will be the 20% of the maximum
    % response any cell in the recording


    for i = 1:length(fieldnames(CALCIUM.BASELINE_TO_PEAK))
        temp5 = fieldnames(CALCIUM.BASELINE_TO_PEAK);
        temp6 = char(temp5(i));
        if isempty(CALCIUM.BASELINE_TO_PEAK.(temp6))
        else
            maxcalcium_response(i) = max(CALCIUM.BASELINE_TO_PEAK.(temp6));
        end
    end

    MAXCALCIUM_RESPONSE = max(maxcalcium_response);
    CALCIUM.MAXCALCIUM_RESPONSE = MAXCALCIUM_RESPONSE;
    CALCIUMimg('save',CALCIUM,[],list,nfish);

    clear MAXCALCIUM_RESPONSE maxcalcium_response temp5 temp6
end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                          CALCIUM IMAGING + SWIMMING(DEEPLABCUT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---- RESPONSE STRUCTURE (WHICH CELL IS RESPONDING TO EACH SWIMMING EPISODE AND WICH CALCIUM TRANSIENT IS)
    %------------------------- WINDOWING (USING X SECONDS BEFORE THE PEAK)

clear 
user_settings

for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;
% ----------------------------------------------------------------------------

    % Because the Calcium signal is slower, the time of the peak of the calcium response
    % may not be inside the time of the peaks of the swimming episodes, we
    % have to add a constant in seconds just to be safe

    clear WINDOW_X WINDOW_X_SWIM WINDOW_Y WINDOW_Y_SWIM

    % We have a period of x seconds after the finish of the swimming episode for the
    % calcium signal to appear
    cte = 3; % seconds.

    % Having the CALCIUM peak we will take a retrospectic window period of x secs
    restric = 3; % seconds.

    % You can used the ONSET of the response or the PEAKTIME
    % TEMP = CALCIUM.ONSET_X;
    TEMP = CALCIUM.AMPLITUDE_X;

    %threshold_calcium_response = 0.10;

    for i = 1:length(fieldnames(CALCIUM.AMPLITUDE_X)) % Loop across all the cells (number of cells)

        yyaxis right
        plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,i)-CALCIUM.baseline_amplitude(i),'b'),
         yyaxis left
            plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x-CALCIUM.baseline_swimming,'g','LineWidth',1);

           

        temp1 = fieldnames(CALCIUM.AMPLITUDE_X);
        myfield1 = char(temp1(i)); % extracting the fieldnames of each cell of the loop
        RESPONSE.(myfield1) = zeros(1,length(fieldnames(CALCIUM.EPISODE_X)));  % prelocating wih the number of swim episodes,
        % if a cell is not responding, it will have a 0.

        for j=1:length(CALCIUM.AMPLITUDE_X.(myfield1)) % number of time point responses of one cell (ONSET)

            % if CALCIUM.BASELINE_TO_PEAK.(myfield1)(j)>threshold_calcium_response*CALCIUM.MAXCALCIUM_RESPONSE


            leftrstr = TEMP.(myfield1)(j)-restric ;  % restriction left
                window = [leftrstr TEMP.(myfield1)(j)]; % this is the retrospectic window that we will
                % use to look at the swimming

                for h= 1:length(fieldnames(CALCIUM.EPISODE_X)) % Loop across the number of episodes
                    temp2 = fieldnames(CALCIUM.EPISODE_X);
                    myfield2 = char(temp2(h)); % extracting the name of each episode

                    unico = CALCIUM.EPISODE_X.(myfield2).*sr; %in seconds
                    ctevalues = unico(end) + cte; %we add cts number of seconds
                    safetime = [unico(1) ctevalues]; % the new vector time just to see in
                    % calcium wave is inside the period of the swimming plus the
                    % cts seconds


                    % THIS is for the RESPONSE STRUCTURE

                    if TEMP.(myfield1)(j)>=safetime(1) && TEMP.(myfield1)(j)<=safetime(2)
                        % If is 1 means that the code founded the location of one
                        % peak of the cell in the period created, so this means the
                        % cell is active at the same time (+ cts seconds) that
                        % swimming is ocurring
                       if RESPONSE.(myfield1)(h)~=0
                       else 
                        RESPONSE.(myfield1)(h) = j;  %  we create this structura that will tell us
                        % how many times a swimming precedes a peak. For
                        % example if two swimming episodes preceded the
                        % third peak it will be [3,3]
                       end 


                        % Here can be added that if RESPONSE.(myfield1)(h)==1 then
                        % i+1 or something like that because maybe there is more
                        % than one response to the same swimming episode.

                        % THIS is for the WINDOW
                        if window(1)>safetime(end) || window(end)<safetime(1) % when there is not overlap btw
                            % the vectors means that there is not timepoints overlaping
                            % so the previous to the calcium activity there is not
                            % swimming response
                        else
                            w1 = CALCIUM.EPISODE_X.(myfield2).*sr>window(1); % elements in the swimming episode bigger
                            % then the window 1
                            w2 = CALCIUM.EPISODE_X.(myfield2).*sr<window(2); % elements in the swimming episode smaller
                            % then the window 2
                            WINDOW_X.(myfield1).(myfield2) = (CALCIUM.EPISODE_X.(myfield2).*sr).*w1.*w2;  % this are only the timepoints of the swimming included in the window
                            % this are only the VALUES of the swimming included in the window
                            WINDOW_Y.(myfield1).(myfield2) = (CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield2)).*w1.*w2;
                            clear w1 w2

                            %%%%Ploting for visualization
                            yyaxis left
                            hold on
                            safetimevector = linspace(safetime(1),safetime(2),100);
                            A = CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield2)(1).*ones(1,length(safetimevector));
                            scatter(safetimevector,A,'r')
                            if WINDOW_X.(myfield1).(myfield2)==0
                            else
                            xline(nonzeros(WINDOW_X.(myfield1).(myfield2)),'r')
                            end
                             yyaxis right
                             hold on
                             window_vector = linspace(window(1),window(2),100);

                            C = CALCIUM.BASELINE_TO_PEAK.(myfield1)(j).*ones(1,length(window_vector));
                            scatter(window_vector,C,'b')


                            %%%%%%%

                            if isempty(CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield2)) % The same for the swimming
                            else
                                w1 = CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield2)>window(1);
                                w2 = CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield2)<window(2);
                                WINDOW_X_SWIM.(myfield1).(myfield2) = CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield2).*w1.*w2;
                                WINDOW_Y_SWIM.(myfield1).(myfield2) =CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2).*w1.*w2;
                                temp10 = nonzeros(WINDOW_X_SWIM.(myfield1).(myfield2)); 
                                if length(temp10)==1
                                    WINDOW_Y_SWIM_FREQ.(myfield1).(myfield2)=[];
                                else
                                for jjj = 1:length(temp10)-1
                                    WINDOW_Y_SWIM_FREQ.(myfield1).(myfield2)(jjj) =  1/(temp10(jjj+1)... 
                                        - temp10(jjj))*2;


                                end 
                                clear w1 w2
                                end 
                            end


                        end



                    else
                        if RESPONSE.(myfield1)(h)~=0 % becasue the construction of the code is overwritting
                            % the previus iteration, we have to keep the ones
                            % previousl created, this is, if the cell response to
                            % the previous episode of swimming.
                        else
                            RESPONSE.(myfield1)(h) = 0;
                        end





                    end
                end
                %%%
%             end
        end


            saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
                char(CALCIUM.list(nfish,1)),strcat('Windows',myfield1, char(CALCIUM.list(nfish,1)))))
close all
    end

    CALCIUM.WINDOW_X = WINDOW_X;
    CALCIUM.WINDOW_Y = WINDOW_Y;
    CALCIUM. WINDOW_X_SWIM =  WINDOW_X_SWIM;
    CALCIUM. WINDOW_Y_SWIM= WINDOW_Y_SWIM;
    CALCIUM.WINDOW_X_SWIM_FREQ = WINDOW_Y_SWIM_FREQ; 
    CALCIUM.RESPONSE = RESPONSE;

    
    CALCIUMimg('save',CALCIUM,[],list,nfish);

    clearvars -except iiii
    user_settings

end


% TO SUM UP: this code is, taking one cell, taking each of the swimming
% episodes (taking the time when this this happening). Then we will take
% one this times and we will check if is inside the times where any episode
% of swimming is ocurring. if is inside it will put 1 otherwise 0. then it
% will take the next timepoint of the calcium response and it will check if
% is inside any episode of swimming. to not overwrite the previous vector
% created we will keep the ones. at the end, for each cell you will have a
%vector with length equal to the number of swimming episodes, and in each
%index will be a 1 if that cell responded to that specific swimming
%episode.



%% CREATING DATA SHEET

%The first row will be the name of the cell
writecell(name_perloc',fullfile(path.data...
    ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),'Sheet',1,'Range','A1')
% The second row will be the position normalize
writematrix(CALCIUM.Normalized_x_distance,fullfile(path.data...
    ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),'Sheet',1,'Range','A2')

sum1 = 0;
for i = index_order
    temp1 = fieldnames(WINDOW_Y);
    myfield1 = char(temp1(i)); 
    sum1 = length(i) + sum1; % for natural counting (1,2,3,...)
    longitud = 0;
    sum = 0; 
    for h=1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
        temp2 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
        myfield2 = char(temp2(h)); % extracting the name of each episode


        if isfield(WINDOW_Y.(myfield1),(myfield2))
            T([1+sum:length(WINDOW_Y.(myfield1).(myfield2))+sum],1) =...
                WINDOW_Y.(myfield1).(myfield2)
            sum = sum + length(WINDOW_Y.(myfield1).(myfield2))
        end
    end
    T = T(find(T)); % Removinf the zeros
    writematrix(T,fullfile(path.data...
        ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),'Sheet',1,'Range',...
        strcat(char('A'+(sum1-1)),'3'))
    clear longitud sum T

end


% Connect to Excel
Excel = actxserver('excel.application');
% Get Workbook object
WB = Excel.Workbooks.Open(fullfile(path.data...
    ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),0,false);

% Set the color of cell "A1" of Sheet 1 to RED
WB.Worksheets.Item(1).Range(strcat('A2:',strcat(char('A'+(length(name_perloc)-1)),'2'))).Interior.ColorIndex = 3;
% Save Workbook
WB.Save();
% Close Workbook
WB.Close();
% Quit Excel
Excel.Quit();


%% MAXIMUM CALCIUM RESPONSE IN ALL THE RECORDINGS (FOR ALL THE RECORDINGS OF THE SAME ANIMAL)
user_settings
for iiii=1:length(list)
     nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;


MAXCALCIUM_RESPONSE_TOTAL(iiii) = CALCIUM.MAXCALCIUM_RESPONSE; 
end



for iiii=1:length(list)
     nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;

    CALCIUM.MAXCALCIUM_RESPONSE_TOTAL = max(MAXCALCIUM_RESPONSE_TOTAL);


CALCIUMimg('save',CALCIUM,[],list,nfish);
end 









