%% David Madrid . Last Rev 5/04/2022

%%                                 CLEAR

clear
close all
user_settings;
list

%%                              LOADING FISH IN A LOOP

for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                                    SWIMMING (DEEPLABCUT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%                 CUTTING SWIMMING AND CALCIUM WAVES

    % CUT THE SWIMMING TRACE AND CALCIUM WAVE TO REMOVE SEGMENTS WHEN THE
    % ACTIVITY IS NOT CLEAR OR NOT ANALYZABLE

    % Plot for visual checking
    plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x)
    pause(0.5)

    prompt = {'Do you want to cut the trace?'};
    dlgtitle = 'Input';
    dims = [1 35];

    definput ={'N'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    answer = char(answer);

    if isequal(answer,'Y')

        disp(['Click  once at the begining of the segmenet that you want to cut and' ...
            'twice at the end'])

        [x_time,~] = getpts;

        % Set initial and terminal points (x-values). The number of decimals that
        % you see on the plot is rounded. Select that number with that number of decimals.
        XCUT1 = x_time(1); % Initial x value where you want to start cutting
        XCUT2 = x_time(2); % Final x value where you want to stop cutting


        substraction_swim = abs(CALCIUMroiTS.deeplabcut.lasertimesync-XCUT1);
        substraction_CAL = abs(CALCIUMroiTS.diff_perc03.times-XCUT1);
        xi_swim_cutted = find(substraction_swim==min(substraction_swim));
        xi_calcium_cutted = find(substraction_CAL==min(substraction_CAL));
        clear substraction_CAL substraction_swim % Erase temporal variable

        %Next we perform the same for the XCUT2
        %For the swimming
        substraction_swim = abs(CALCIUMroiTS.deeplabcut.lasertimesync-XCUT2);
        substraction_CAL = abs(CALCIUMroiTS.diff_perc03.times-XCUT2);
        xf_swim_cutted = find(substraction_swim==min(substraction_swim));
        xf_calcium_cutted = find(substraction_CAL==min(substraction_CAL));
        clear substraction_CAL substraction_swim % Erase temporal variable


        % At the end of the computations you end with 4 variables (xi_swim_cutted, xf_swim_cutted,
        % xi_calcium_cutted,xf_calcium_cutted) and the cutting points XCUT1 and
        % XCUT2


        % NEXT, we create the new waves, both for the CALCIUM (time vector and data vector)
        % and for the swimming (time vector and data vector). We will overwritted
        % the actual variables, for getting the variables back just rerun the first
        % part of the code and they will be loaded again

        % For the swimming

        CALCIUMroiTS.deeplabcut.lasertimesync = ...
            CALCIUMroiTS.deeplabcut.lasertimesync(xi_swim_cutted:xf_swim_cutted);

        CALCIUMroiTS.deeplabcut.tail3.x = ...
            CALCIUMroiTS.deeplabcut.tail3.x(xi_swim_cutted:xf_swim_cutted);

        % For the CALCIUM

        CALCIUMroiTS.diff_perc03.times = ...
            CALCIUMroiTS.diff_perc03.times(xi_calcium_cutted:xf_calcium_cutted);
        CALCIUMroiTS.diff_perc03.data=...
            CALCIUMroiTS.diff_perc03.data(xi_calcium_cutted:xf_calcium_cutted,:);

        % --------------- Test Ploting ----------------

        plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x)
        hold on
        plot(CALCIUMroiTS.diff_perc03.times,CALCIUMroiTS.diff_perc03.data(:,1))

        % Restarting the time vector (both swimming and CALCIUM)

        % The time vectors will start for their respective sampling rates, so both
        % of them are restarted

        CALCIUMroiTS.deeplabcut.lasertimesync = CALCIUMroiTS.deeplabcut.lasertimesync...
            - CALCIUMroiTS.deeplabcut.lasertimesync(1)+sr;

        CALCIUMroiTS.diff_perc03.times = CALCIUMroiTS.diff_perc03.times...
            -CALCIUMroiTS.diff_perc03.times(1)+CALCIUM.srate;

        % --------------- Test Ploting -------------------------------

        plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x)
        hold on
        plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,4))


        % Saving the cutting points

        CALCIUM.XCUT1 = XCUT1;
        CALCIUM.XCUT2 = XCUT2;

        % Close figures

        close all

    end

    close all

    CALCIUMimg('save',CALCIUM,[],list,nfish);
    CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);

    
    %%              PEAKS IDENTIFICATION AND CLEANING OF THE SWIMMING TRACE


    % ------------ VISUAL CHECKING FO THE PEAKS (SET UP PEAKS THERSHOLD)

    sidepeaks =10; % Initial value of the threshold for the peaks

    % Plot for visual checking
    figure('units','normalized','outerposition',[0 0 1 1])
    findpeaks(CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidepeaks);

    pause(2) % Time to plot the figure

    % This lil prompt will ask if you are satisfy or not, if the str=N,
    % means that you dont want to modify  the threeshold
    prompt = 'PEAKS: If you want to modify the threshold for the peaks type "Y", if you are satisfy type "N" : ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
        close
    else

    end


    % If the answer is str=Y you will enter in a while loop that will allow you
    % to change the value of the threshold until you are satisfy

    while isequal(str,'Y')

        pause(1)

        prompt = {'Enter Threshold:'};
        dlgtitle = 'Input';
        dims = [1 35];

        definput ={sprintf('%.0f',sidepeaks)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        sidepeaks = str2double(answer(1));
        cla
        findpeaks(CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidepeaks);

        pause(2)

        prompt = 'PEAKS: If you want to modify the threshold for the peaks type "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');

        if isempty(str)
            str = 'N';
            close
        end

    end

    % ------------ VISUAL CHECKING FOR THE TROUGHS (SET UP TROUGHS THERSHOLD)

    sidetroughts = 20; % Initial value of the threshold for the troughs

    % Plot for visual checking
    figure('units','normalized','outerposition',[0 0 1 1])
    findpeaks(-CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidetroughts);

    pause(2) % Time to plot the figure

    % This lil prompt will ask if you are satisfy or not, if the str=N,
    % means that you dont want to modify  the threeshold
    prompt = 'TROUGHS: If you want to modify the threshold for the peaks type "Y", if you are satisfy type "N" : ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
        close
    else

    end


    % If the answer is str=Y you will enter in a while loop that will allow you
    % to change the value of the threshold until you are satisfy

    while isequal(str,'Y')

        pause(1)

        prompt = {'Enter Threshold:'};
        dlgtitle = 'Input';
        dims = [1 35];

        definput ={sprintf('%.0f',sidetroughts)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        sidetroughts = str2double(answer(1));
        cla
        findpeaks(-CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidetroughts);
        pause(2)

        prompt = 'TROUGHS: If you want to modify the threshold for the peaks type "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');

        if isempty(str)
            str = 'N';
            close
        end

    end

    close all

    % -------------------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------------------


    %                           FINDING PEAKS AND TROUGHS



    %----------------------    Finding and removing peaks

    figure('units','normalized','outerposition',[0 0 1 1])
    findpeaks(CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidepeaks); hold on
    [pks1,locs1] = findpeaks(CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidepeaks);

    prompt = 'PEAKS: If you want to remove some peaks "Y", if you are satisfy type "N" : ';

    str = input(prompt,'s');
    if isequal(str,'N')
        str = 'N';
        close

    elseif isequal(str,'Y')

        disp('Clik on the peaks that you want to remove, double click in the last one.Push delete to remove a point')
        [x_time,~] = getpts;


        % Finding closest match
        removing = NaN(1,length(x_time)); % Prelocate
        for i=1:length(x_time)
            removing(i) = find(abs(locs1-x_time(i))==min(abs(locs1-x_time(i))));
        end

        % Removing this index values from the peaks and locs vectors
        pks1(removing) = [];
        locs1(removing) = [];

        % Testing

        scatter(locs1,pks1,30,'filled','g')
        pause(0.5)

        prompt = 'Now you can zoom in. If you want to keep removing points type "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');

        if isempty(str)
            str = 'N';
            close
        end


        while isequal(str,'Y')
            [x_time,~] = getpts;
            % Finding closest match
            removing = NaN(1,length(x_time)); % Prelocate
            for i=1:length(x_time)
                removing(i) = find(abs(locs1-x_time(i))==min(abs(locs1-x_time(i))));
            end

            % Removing this index values from the peaks and locs vectors
            pks1(removing) = [];
            locs1(removing) = [];
            cla
            findpeaks(CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidepeaks); hold on
            scatter(locs1,pks1,30,'filled','g')
            pause(0.5)
            prompt = 'Now you can zoom in. If you want to keep removing points type "Y", if you are satisfy type "N" : ';
            str = input(prompt,'s');

            if isempty(str)
                str = 'N';
                close
            end
        end
    end


    clear x_time y_time removing
    close all


    % -------------     Finding  and removing troughs


    figure('units','normalized','outerposition',[0 0 1 1])
    findpeaks(-CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidetroughts); hold on
    [pks2,locs2] = findpeaks(-CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidetroughts);

    prompt = 'TROUGHS: If you want to remove some peaks "Y", if you are satisfy type "N" : ';
    str = input(prompt,'s');

    if isequal(str,'N')
        str = 'N';
        close

    elseif isequal(str,'Y')
        disp('Clik on the troughs that you want to remove, double click in the last one. Push delete to remove a point')
        [x_time,~] = getpts;


        % Finding closest match
        removing = NaN(1,length(x_time)); % Prelocate
        for i=1:length(x_time)
            removing(i) = find(abs(locs2-x_time(i))==min(abs(locs2-x_time(i))));
        end
        % Removing this index values from the peaks and locs vectors
        pks2(removing) = []; locs2(removing) = [];

        % Testing

        scatter(locs2,pks2,30,'filled','g')
        pause(0.5)

        prompt = 'If you want to keep removing points type "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');

        if isempty(str)
            str = 'N'; % you are done
            close
        end


        while isequal(str,'Y') % if the answe is yes you enter in a while loop to elimite other peaks
            [x_time,~] = getpts;
            % Finding closest match
            removing = NaN(1,length(x_time)); % Prelocate
            for i=1:length(x_time)
                removing(i) = find(abs(locs2-x_time(i))==min(abs(locs2-x_time(i))));
            end

            % Removing this index values from the peaks and locs vectors
            pks2(removing) = [];
            locs2(removing) = [];
            cla
            findpeaks(-CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidetroughts); hold on
            scatter(locs2,pks2,30,'filled','g')
            pause(0.5)

            prompt = 'If you want to keep removing points type "Y", if you are satisfy type "N" : ';
            str = input(prompt,'s');

            if isempty(str)
                str = 'N'; % you are done
                close
            end
        end
    end


    clear x_time y_time removing
    close all


    % Ploting everything together

    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(211)
    title('Peaks')
    findpeaks(CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidepeaks); hold on
    scatter(locs1,pks1,30,'filled','g')

    subplot(212)
    title('Troughs')
    findpeaks(-CALCIUMroiTS.deeplabcut.tail3.x,'MinPeakProminence',sidetroughts); hold on
    scatter(locs2,pks2,30,'filled','g')
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Accurate_peaks_troughs',char(CALCIUM.list(nfish,1)),'.fig')))
    pause(2)
    close all
    clear prompt str


    % Saving the values of the threshold for peaks and troughs and values

    CALCIUM.sidepeaks_swim = sidepeaks;
    CALCIUM.sidetroughts_swim  = sidetroughts;
    CALCIUM.pks1_swim  = pks1;
    CALCIUM.pks2_swim  = pks2;
    CALCIUM.locs1_swim  = locs1;
    CALCIUM.locs2_swim  = locs2;

    CALCIUMimg('save',CALCIUM,[],list,nfish);

    %%                        MANTEINING THE SING OF PEAKS AND TROUHTS

    time_vector_sing = sort ( [CALCIUM.locs1_swim ; CALCIUM.locs2_swim] );
    for i=1:length(time_vector_sing)
        if any(time_vector_sing(i) == CALCIUM.locs1_swim)
            temp_mantain =  any(time_vector_sing(i) == CALCIUM.locs1_swim); 
            sign_retention (i) = CALCIUM.pks1_swim(temp_mantain);
             

        elseif any(time_vector_sing(i) == CALCIUM.locs2_swim) 
            temp_mantain =  any(time_vector_sing(i) == CALCIUM.locs2_swim);
            sign_retention (i) = CALCIUM.pks2_swim(temp_mantain); 

    end 
    end 
    
    % We will take just the fact that is positive (1) or negative (-1)

  for i=1:length(sign_retention )
      if sign_retention (i)>0
          sign_retention (i) = 1;
      else 
          sign_retention (i) = -1;
      end 
  end 

   CALCIUM.sign_retention = sign_retention;
   CALCIUMimg('save',CALCIUM,[],list,nfish);

    %%                      AISLATION OF SWIMMING EPISODES AND BASELINE PER EACH EPISODE (MANUAL)

    % Select manually a segment that contains each episode. So each two consectuve number in
    % vector created will correspondto the start and end of the segmenet and
    % the swimming episode will be contain in that segment.
    figure('units','normalized','outerposition',[0 0 1 1])
    
    plot(CALCIUMroiTS.deeplabcut.tail3.x), hold on
    title('EXTRACTING SWIMMING EPISODES')

    disp('You can zoom in now')

    prompt = 'Do you want to select another episodes "Y", if you are satisfy type "N" : ';
    str = input(prompt,'s');
    sum = 0;

    while isequal(str,'Y')

        disp('click a two times. Inside two points you must find a swimming episode.')
        [x_time,~] = getpts; %@

        points(1,[1+sum 2+sum]) = x_time; 
        prompt = 'Do you want to select another episodes "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');
        clear x_time
        sum = sum + 2;

    end
    
   

%     FOR THE BASELINE OF EACH EPISODE

 
   
    plot(CALCIUMroiTS.deeplabcut.tail3.x), hold on
     title('EXTRACTING BASELINE OF EACH EPISODES')

    disp('Select a small segment right before the first peak of each episode')
    disp('You can zoom in now')

    prompt = 'Do you want to select another baseline "Y", if you are satisfy type "N" : ';
    str = input(prompt,'s');
    sum = 0;

    while isequal(str,'Y')

        disp('click a two times. Inside two points you must find a swimming episode.')
        [x_time_B,~] = getpts; %@

        points_B(1,[1+sum 2+sum]) = x_time_B; 
        prompt = 'Do you want to select another baseline "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');
        clear x_time_B
        sum = sum + 2;

    end

   


    x_time= sort(points); % Swimming episodes
    x_time_B = sort(points_B); % baselines for each episodes


    % We need to have even number in the x_time variable: Control test
    if rem(length(x_time), 2) == 0 && rem(length(x_time_B), 2)==0 &&...
            length(x_time)==length(x_time_B)

    else
        error('The number of points is not even or more beaselines than episodes or the opposite')
    end



 % Finfing the indexs that correspond with the points we just got

    for i=1:length(x_time_B)
        substraction = abs(CALCIUMroiTS.deeplabcut.lasertimesync-x_time_B(i)*sr);
        xi_baseline_swim (i) = find(substraction==min(substraction));
        clear substraction % Erase temporal variable
    end

% Computating the avrg
 sum = 0; 
for i=1:length(x_time_B)/2
    baseline_swimming (i)=...
        mean(CALCIUMroiTS.deeplabcut.tail3.x([xi_baseline_swim(1+sum):xi_baseline_swim(2+sum)])); %avrg

    std_baseline (i) = std(CALCIUMroiTS.deeplabcut.tail3.x([xi_baseline_swim(1+sum):xi_baseline_swim(2+sum)])); %std
   sum = sum + 2;
end 
    % Next we will select the locs1 and locs2 that are within each segment and
    % we sort the elements


sum = 0; 
sum1 = 0; 
    for k=1:2:length(x_time)
        my_field1 = strcat('swm_',num2str(k)); % Creating dinamicly the
        % name of the field
        my_field2 = strcat('swm_',num2str(k));
        index1 = find((CALCIUM.locs1_swim>x_time(k) & CALCIUM.locs1_swim<x_time(k+1))); % This will
        % find  all the values of CALCIUM.locs1_swim inside each interval
        index2 = find((CALCIUM.locs2_swim>x_time(k) & CALCIUM.locs2_swim<x_time(k+1))); % This will
        % find  all the values of CALCIUM.locs2_swim inside each interval
        EPISODE_X.(my_field1) = sort([CALCIUM.locs1_swim(index1);CALCIUM.locs2_swim(index2)]); % we put
        %everything in a column and sort it.
        EPISODE_Y.(my_field2) = CALCIUMroiTS.deeplabcut.tail3.x(EPISODE_X.(my_field1));
        BASELINE_PER_EPISODE.(my_field2) =  baseline_swimming(k-sum); 
        temp_var = sort([CALCIUM.locs1_swim(index1);CALCIUM.locs2_swim(index2)]);
        temp_var_2 = linspace(xi_baseline_swim(1 + sum1), xi_baseline_swim(2 + sum1),100); 
        plot(temp_var_2,baseline_swimming(k-sum)*ones(1,length(temp_var_2)),'R','LineWidth',3); 
        plot(temp_var,baseline_swimming(k-sum).*ones(size(temp_var)),'g','LineWidth',6);
        
        clear temp_var temp_var_2
        sum = sum + length(sum); 
        sum1 = sum1 +2; 
    end

    % Now the variable EPISODES has a number of variables that correspond with
    % each swimming episode


    %Saving

    CALCIUM.EPISODE_X = EPISODE_X;
    CALCIUM.EPISODE_Y = EPISODE_Y;
    CALCIUM.BASELINE_PER_EPISODE = BASELINE_PER_EPISODE;
    CALCIUMimg('save',CALCIUM,[],list,nfish);
    pause(2)
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Baseline_swim_per_episode',char(CALCIUM.list(nfish,1)),'.fig')))

    close all
    clear xi_baseline_swim baseline_swimming  std_baseline

    %%
    % EXTRACTING THE GENERAL BASELINE FOR SWIMMING
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x)
    title('EXTRACTING GENERAL BASELINE')
    disp(['Click once at the begining of the segment that you want to use ' ...
        'as a baseline, and twice at the end of it.'])
    [x_time,~] = getpts; %@
    cla

    % Control just to be sure we select 2 points
    if length(x_time)>2
        error('You click more the two points')
    elseif length(x_time)<2
        error('You click just the first point')
    end

    % Finfing the indexs that correspond with the points we just got
    for i=1:size(x_time,1)
        substraction = abs(CALCIUMroiTS.deeplabcut.lasertimesync-x_time(i,1));
        xi_baseline_swim (i) = find(substraction==min(substraction));
        clear substraction % Erase temporal variable
    end

    % Computating the avrg
    baseline_swimming=...
        mean(CALCIUMroiTS.deeplabcut.tail3.x([xi_baseline_swim(1):xi_baseline_swim(2)])); %avrg

    std_baseline = std(CALCIUMroiTS.deeplabcut.tail3.x([xi_baseline_swim(1):xi_baseline_swim(2)])); %std

    %Visual test
    plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x); hold on
    plot(CALCIUMroiTS.deeplabcut.lasertimesync,baseline_swimming*ones(size(CALCIUMroiTS.deeplabcut.lasertimesync)),'g','LineWidth',4);
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('General_Baseline_swim',char(CALCIUM.list(nfish,1)),'.fig')))

    % Saving
    CALCIUM.baseline_swimming = baseline_swimming;
    CALCIUM.std_baseline_swimming = std_baseline; 
    CALCIUMimg('save',CALCIUM,[],list,nfish);
    pause(2)
    close all

    %%              SWIMMING POWER (USING THE AMPLITUDE TO THE BASELINE FOR EACH EPISODE)

    % We will use the baseline per EACH EPIOSDE. If you want to use the
    % general baseline like tempor

      tempor = CALCIUM.BASELINE_PER_EPISODE;
%     tempor = CALCIUM.baseline_swimming;

    for i=1:length(fieldnames(CALCIUM.EPISODE_Y)) %Loop across all the swimming episodes
        temp = fieldnames(CALCIUM.EPISODE_Y); % We use the "Y" because we need the value of the peaks
        %for each episode
        my_field = char(temp(i)); % In each loop this variables will take the name of
        % one of the fields (one of the swimming episodes)
        EPISODE_Y_BASELINE_AMP.(my_field) = abs(CALCIUM.EPISODE_Y.(my_field)-tempor.(my_field)); % we create
        % this new structure that will have positive amplitudes from the
        % baseline to the peaks
        EPISODE_Y_BASELINE_AMP_AVGR.(my_field) = mean(EPISODE_Y_BASELINE_AMP.(my_field));% we create
        % this new structure that will have averages from the
        % baseline to the peaks

        %Normalize the the amplitude to the baseline
        maxim(i) = max(EPISODE_Y_BASELINE_AMP.(my_field)); % This variable contain all the max for each
        %swimming episode
        maxmaxim = max(maxim); % this variable contain the maximal amplitude during the whole recording
        EPISODE_Y_BASELINE_AMP_NORM.(my_field) = EPISODE_Y_BASELINE_AMP.(my_field)./maxmaxim; %This struct
        % is the same that EPISODE_Y_BASELINE_AMP.(my_field) but normalized.
        EPISODE_Y_BASELINE_AMP_AVGR_NORM.(my_field) = mean(EPISODE_Y_BASELINE_AMP_NORM.(my_field)); %This struct
        % is the same that EPISODE_Y_BASELINE_AMP_AVGR.(my_field) but normalized.

    end

    % Saving
    CALCIUM.EPISODE_Y_BASELINE_AMP = EPISODE_Y_BASELINE_AMP;
    CALCIUM.EPISODE_Y_BASELINE_AMP_AVGR = EPISODE_Y_BASELINE_AMP_AVGR;
    CALCIUM.EPISODE_Y_BASELINE_AMP_NORM = EPISODE_Y_BASELINE_AMP_NORM;
    CALCIUM.EPISODE_Y_BASELINE_AMP_AVGR_NORM = EPISODE_Y_BASELINE_AMP_AVGR_NORM;

    CALCIUMimg('save',CALCIUM,[],list,nfish);
    clear tempor


    %%                              EPISODES WITH SIGN
    
    acum = 0; 
    for i = 1:length(fieldnames(CALCIUM.EPISODE_X))
        myfield = fieldnames(CALCIUM.EPISODE_X); 
        temp_sing = char(myfield(i)); 

        numer = length(CALCIUM.EPISODE_Y_BASELINE_AMP.(temp_sing)); 
        
        CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(temp_sing) = CALCIUM.sign_retention ([acum+1:acum + numer]); 
        acum = numer + acum;
    end

    CALCIUMimg('save',CALCIUM,[],list,nfish);

    %%                              REAL SWIMMING ANALYSIS (ALTERNATION)



    % We will use the baseline per EACH EPIOSDE. If you want to use the
    % general baseline like tempor

      tempor = CALCIUM.BASELINE_PER_EPISODE;
%     tempor = CALCIUM.baseline_swimming;

    for k=1:length(fieldnames(CALCIUM.EPISODE_X)) % we will do it across swimming episodes
        temp1 = fieldnames(CALCIUM.EPISODE_X);
        myfield1 = char(temp1(k));

        % REMOVING POINTS FURTHER AWAY BASED ON SMALL FREQUENCY

        timepoints = CALCIUM.EPISODE_X.(myfield1);

        %Forward identification
        for i=1:length(timepoints)-1
            distance_frames(i) = timepoints (i+1)-timepoints (i);
            freq(i) = 1/(distance_frames(i)*(sr));
        end

        threshold = 4; %Remove date points that are really far away (small frequency), posible
        %different episodes f swimming

        v=zeros(1,length(timepoints));

        for i=1:length(timepoints)-1
            if freq(i)>threshold
                v (i) = 1;
                v (i+1) = 1;

            else
                v(i+1) = 0;
            end
        end
        clear freq distance_frames threshold

        v1 = find(v);
        % THIS VECOR (V) HAVE THE LOGICAL VALUE OF THE TIMEPOINTS WHERE THE FREQUENCY BTW
        % HALF OF A CYCLE IS BIGGER THAN A THRESHOLD

        if isempty(v1) % If there is no swimming there will be no index and all the rest of the code will give
            % us error, so if v1 is empty means that in this episode there
            % is not swimming and we can jumo to the next one


            CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield1) =  [];
            CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield1) = [];
             CALCIUM.FREQ_SWIM_X.(myfield1) = [];
                CALCIUM.FREQ_SWIM_Y.(myfield1) = [];
            CALCIUMimg('save',CALCIUM,[],list,nfish);

        else % if is not empty mean there is swimming and we can continue with the code


            plot(CALCIUM.EPISODE_X.(myfield1)*sr,tempor.(myfield1)*ones(size(CALCIUM.EPISODE_X.(myfield1))),'g','LineWidth',4);hold on
            plot(CALCIUM.EPISODE_X.(myfield1)*sr,CALCIUM.EPISODE_Y.(myfield1),'LineWidth',3);
            xline(timepoints(v1)*(sr),'LineWidth',1)

            % REMOVE CONSECUTIVE PEAKS OR TROUGTS WITH SMALL DIFERENCES

            % % Equal number of data points in peaks, if is not equal this will add zeros
            % % at the end of the shorter vector
            %
            % if length(CALCIUM.pks1_swim)==length(CALCIUM.pks2_swim)
            %     return
            % elseif length(CALCIUM.pks1_swim)>length(CALCIUM.pks2_swim)
            %     zeropading = zeros(1,length(CALCIUM.pks1_swim)-length(CALCIUM.pks2_swim));
            %     CALCIUM.pks2_swim = [CALCIUM.pks2_swim; zeropading'];
            % elseif length(CALCIUM.pks1_swim)<length(CALCIUM.pks2_swim)
            %     zeropading = zeros(1,length(CALCIUM.pks2_swim)-length(CALCIUM.pks1_swim));
            %     CALCIUM.pks1_swim= [CALCIUM.pks1_swim; zeropading'];
            % end

            % Next we have to organize time-wise the peaks and troughs and puting them
            % in orde on a vector. For that i use the conditionals to first detect if
            % the time of and first event (peak or trougts) is found in the vector
            % location of the peaks (locs1) or of the trougts(locs2). them i will take
            % the peak or the trougt corresponding to that time point and i will place
            % it in a new vector. Them i will do it with the second timepoint and so on
            % Now the vector event_vector contains the values of the peaks and trougts
            % cronologally order


            Amplitud_threshold = 40;

            % for i = 1:length(timepoints)
            %     if  ismember(timepoints(i),CALCIUM.EPISODE_X.(myfield1))
            %         event_vector(i) =CALCIUM.EPISODE_Y.(myfield1)(find(timepoints(i)==CALCIUM.EPISODE_X.(myfield1)));
            %     elseif ismember(timepoints(i),CALCIUM.locs2_swim)
            %         event_vector(i) =CALCIUM.pks2_swim(find(timepoints(i)==CALCIUM.locs2_swim));
            %     end
            % end


            % Remove elements with small diferences (less than amplitude_threshold)

            remove_amplitude=NaN(1,length(timepoints));
            for i=1:length(timepoints)-1
                if abs(abs(CALCIUM.EPISODE_Y.(myfield1)(i+1))-abs(CALCIUM.EPISODE_Y.(myfield1)(i)))>Amplitud_threshold
                    if remove_amplitude (i) == 0
                        remove_amplitude (i+1) = 1;
                    else
                        remove_amplitude (i) = 1;
                        remove_amplitude (i+1) = 1;
                    end
                else
                    remove_amplitude (i+1) = 0;
                end
            end

            v2 = find(remove_amplitude);
            xline(timepoints(v2)*(sr),'r--','LineWidth',1)

            ii = intersect(v1,v2);
            xline(timepoints(ii)*(sr),'g--','LineWidth',1)

            clear Amplitud_threshold  remove_amplitude

            % REMOVE JUST TAILBEATING

            tail_beating=zeros(1,length(timepoints));
            tail_beating(1)=1;
            std_used = 30; % how many standars deviation above and below the baseline

            % Remove the values that are consecutivetly way below the baseline
            baseline = tempor.(myfield1);
            topbaseline = baseline + std_used*std_baseline;
            bottombaseline = baseline - std_used*std_baseline;
            for  i=1:length(timepoints)-1
                if abs(CALCIUM.EPISODE_Y.(myfield1)(i))<bottombaseline && abs(CALCIUM.EPISODE_Y.(myfield1)(i+1))<topbaseline
                    if tail_beating(i)==1
                        tail_beating(i)=1;
                    else
                        tail_beating(i)=0;
                    end

                elseif abs(CALCIUM.EPISODE_Y.(myfield1)(i))>topbaseline && abs(CALCIUM.EPISODE_Y.(myfield1)(i+1))>bottombaseline
                    if tail_beating(i)==1
                        tail_beating(i)=1;
                    else
                        tail_beating(i)=0;
                    end
                else

                    tail_beating(i+1)=1;
                end
            end

            v3= find(tail_beating);
            xline(timepoints(v3)*(sr),'b--','LineWidth',1)


            jj = intersect(ii, v3);
            xline(timepoints(jj)*(sr),'m--','LineWidth',3)

            clear tail_beating

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % CONTINUITY OF THE DATA
            %Remove some datapoints that are far away and solitary

            diff_distance = 0.20;

            datapoints = timepoints(jj)*(sr);
            if length(datapoints)==1 %means there is no swimming
                 CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield1) =  [];
                CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield1) = [];
                CALCIUM.FREQ_SWIM_X.(myfield1) = [];
                CALCIUM.FREQ_SWIM_Y.(myfield1) = [];
                CALCIUMimg('save',CALCIUM,[],list,nfish);
                close all

            else

                for i = 1:length (datapoints)-1
                    diff_datapoints(i) = abs(datapoints(i)-datapoints(i+1));
                    if diff_datapoints(i)<diff_distance
                        v4(i)=1
                        v4(i+1) = 1
                    else
                        v4(i+1)=0

                    end
                end
                v5= find(v4);





                xline(datapoints(v5),'g--','LineWidth',5)

                Final_points_swimming = datapoints(v5);


                saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
                    char(CALCIUM.list(nfish,1)),strcat('Swimming_',myfield1,'_',char(CALCIUM.list(nfish,1)),'.fig')))

                close all

                CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield1) =  Final_points_swimming ;
                CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield1) = abs(CALCIUM.EPISODE_Y.(myfield1)(v5)-tempor.(myfield1));
                CALCIUMimg('save',CALCIUM,[],list,nfish);

                % here you can add CALCIUM.FINAL_POINTS_SWIMMING_Y_BASELINE_AMP ,
                % CALCIUM.FINAL_POINTS_SWIMMING_Y_BASELINE_AMP_AVRG, ...ETC SIEMPRE
                % INDEXANDO CON V5


                clear v v1 v2 ii v3 jj v4 datapoints v5
                close all



                % EXTRACTION OF THE INSTANTANEUS FREQUENCY HALF CYCLE
                for i=1:length(Final_points_swimming)-1
                    distance_frames_sw(i) = Final_points_swimming(i+1)-Final_points_swimming(i);
                    freq_sw(i) = 1./(distance_frames_sw(i))*2;
                end


                % Finding the indexs of the points the frequencies which btw different episodes are higher than 5 hz
                above_5 = find(freq_sw>3); % like this we remove the frq btw swimming episodes that will be lower

                % POTING INSTANTANEUS FREQUENCY VS AMPLITUD CALCIUM IMAGING
                for i=1:length(Final_points_swimming)-1
                    middpoint(i) = (Final_points_swimming(i+1)+Final_points_swimming(i))/2;
                end

                CALCIUM.FREQ_SWIM_X.(myfield1) = middpoint;
                CALCIUM.FREQ_SWIM_Y.(myfield1) = freq_sw(above_5);
                CALCIUMimg('save',CALCIUM,[],list,nfish);
                clear distance_frames_sw freq_sw  above_5 middpoint

            end
        end


    end



%%
    close all
    clearvars -except iiii
    user_settings
    
end

%%%%%% NORMALIZING ALL THE SWIMMING EPISODES AND TAIL BEATING OF THE ANIMAL

% TAIL BEATING

for  i=1:length(list)

    nfish = i;
    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish); 

    for ii = 1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
        temp = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
        my_field = char(temp(ii));

        TB_amplitude(ii) = max(CALCIUM.EPISODE_Y_BASELINE_AMP.(my_field));
    end 
    amplitude_animal_TB(i) = max(TB_amplitude); 
end
    NMAX_TAILBEATING = max( amplitude_animal_TB); 
    

% REAL SWIMMING

for  i=1:length(list)

    nfish = i;
    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish); 

    for ii = 1:length(fieldnames(CALCIUM.FINAL_POINTS_SWIMMING_Y))
        temp = fieldnames(CALCIUM.FINAL_POINTS_SWIMMING_Y);
        my_field = char(temp(ii));
        if isempty(CALCIUM.FINAL_POINTS_SWIMMING_Y.(my_field))
            swimming_amplitude(ii) = 0;
        else 

       swimming_amplitude(ii) = max(CALCIUM.FINAL_POINTS_SWIMMING_Y.(my_field));
        end 
    end 
    amplitude_animal_SW(i) = max(swimming_amplitude); 
end
    NMAX_SWIMMING = max(amplitude_animal_SW); 


    % FREQUENCY 
    for  i=1:length(list)

        nfish = i;
        [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);

        for ii = 1:length(fieldnames(CALCIUM.FREQ_SWIM_Y))
            temp = fieldnames(CALCIUM.FREQ_SWIM_Y);
            my_field = char(temp(ii)); 
            if isempty(CALCIUM.FREQ_SWIM_Y.(my_field))
            swimming_FREQ(ii) = 0;
            else 
                swimming_FREQ(ii) = max(CALCIUM.FREQ_SWIM_Y.(my_field)); 
            end 

        end
       freq_animal_sw(i) =  max(swimming_FREQ);
    end
    NMAX_FREQ = max(freq_animal_sw); 


% SAVING IN ALL CALCIUM FILES THIS VALUES


    for  i=1:length(list)

        nfish = i;
        [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);

        CALCIUM.NMAX_SWIMMING = NMAX_SWIMMING;
        CALCIUM.NMAX_TAILBEATING = NMAX_TAILBEATING; 
        CALCIUM.NMAX_FREQ=NMAX_FREQ; 
        CALCIUMimg('save',CALCIUM,[],list,nfish);
    end

clear all
close all

