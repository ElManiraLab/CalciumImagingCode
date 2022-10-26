%% David Madrid . Last Rev 5/04/2022

%%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
% ---- all episodes inside the retrospective window
clear
user_settings
for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;


    % will be consider peaks of the calcium waves every wave that is at least
    % 20% amplitude of the bigger amplitud in the recording
    sum=0;  % for counting 1,2,3... in every loop independing of the index

    for i = CALCIUM.index_order % neurons order from medial to lateral. This variable is below

        temp1 = fieldnames(CALCIUM.AMPLITUDE_Y);
        myfield1 = char(temp1(i)); % extracting the name of each neuron
        sum = length(i) + sum;

        for h=1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
            temp2 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
            myfield2 = char(temp2(h)); % extracting the name of each episode
            index=CALCIUM.RESPONSE.(myfield1)(h); % This index will tell us if the cell is not responding
            % to the current episode and also at wich CALCIUM transient is
            % couple the response

            if any(index)         %isfield(CALCIUM.WINDOW_Y,(myfield1)) % Alternative for the if

                % ------------- Next, ploting is for the tail beating

                if isfield(CALCIUM.WINDOW_Y.(myfield1),(myfield2)) % if the field exits, plot the values of the amplitudes
                    %of the swimming included in the window of that peak
                    a=find(CALCIUM.WINDOW_Y.(myfield1).(myfield2));
                    c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                    b=(c.*ones(1,length(CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                    subplot(311)

                    % Ploting all the amplitude of each episode

                    scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)))',...
                        CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)./ CALCIUM.NMAX_TAILBEATING,80,b','filled'), hold on


                end
                clear a b c

                % Next is for the swimming amplitude
                if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
                    if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
                        a=find(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2));
                        c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                        b=(c.*ones(1,length(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                        subplot(312)

                        % Ploting all the amplitude of each episode
                        scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)))',...
                            CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)./  CALCIUM.NMAX_SWIMMING,80,b','filled'), hold on


                        % Next is for the swimming freq
                        if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))

                            if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))

                                a=find(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2));
                                if isempty(a)
                                else
                                    b=(c.*ones(1,length(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)(a))));

                                    subplot(313)
                                    % Ploting all the amplitude of each episode
                                    scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)(a)))',...
                                        CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)',80,b','filled'), hold on
                                end

                            else
                            end
                        else
                        end

                    end
                end
            else
            end
        end
    end
end

subplot(311)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Tail beating')

subplot(312)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the swimming normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Swimming')

subplot(313)
xlim([-1.5 1.5])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Freq half cycle','FontSize',13)
title('Frequency')

sgtitle('ALL POINTS FOR EACH EPISODE')

saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
    strcat('Resume_all_points_episode')))
close all




%%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
% ---- max of all episodes inside the retrospective window
clear
user_settings
for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;


    % will be consider peaks of the calcium waves every wave that is at least
    % 20% amplitude of the bigger amplitud in the recording
    sum=0;  % for counting 1,2,3... in every loop independing of the index
    for i = CALCIUM.index_order % neurons order from medial to lateral. This variable is below

        temp1 = fieldnames(CALCIUM.AMPLITUDE_Y);
        myfield1 = char(temp1(i)); % extracting the name of each neuron
        sum = length(i) + sum;

        for h=1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
            temp2 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
            myfield2 = char(temp2(h)); % extracting the name of each episode
            index=CALCIUM.RESPONSE.(myfield1)(h); % This index will tell us if the cell is not responding
            % to the current episode and also at wich CALCIUM transient is
            % couple the response
            if any(index)         %isfield(CALCIUM.WINDOW_Y,(myfield1)) % Alternative for the if

                % Next is fot the tail beating
                if isfield(CALCIUM.WINDOW_Y.(myfield1),(myfield2)) % if the field exits, plot the values of the amplitudes
                    %of the swimming included in the window of that peak
                    a=find(CALCIUM.WINDOW_Y.(myfield1).(myfield2));
                    c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                    b=(c.*ones(1,length(CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                    subplot(311)


                    % PLoting max of each episode:

                    scatter(CALCIUM.Normalized_x_distance(sum),...
                        max(CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a))./ CALCIUM.NMAX_TAILBEATING,80,b(1),'filled'), hold on

                end
                clear a b c

                % ------------- Next, ploting is fot the tail beating
                if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
                    if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
                     
                        a=find(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2));
                        c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                        if isempty(a)
                        else
                        
                        b=(c.*ones(1,length(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                        subplot(312)


                        % PLoting max of each episode:
                        scatter(CALCIUM.Normalized_x_distance(sum),...
                            max(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a))./  CALCIUM.NMAX_SWIMMING,80,b(1),'filled'), hold on
                        end 
                        % Next is for the swimming freq
                        if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))

                            if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))


                                a=find(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2));
                                if isempty(a)
                                else
                                    b=(c.*ones(1,length(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)(a))));

                                    subplot(313)

                                    % PLoting max of each episode:
                                    scatter(CALCIUM.Normalized_x_distance(sum),...
                                        max(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)),80,b(1),'filled'), hold on
                                end

                            else
                            end
                        else
                        end
                    end
                end
            else
            end
        end
    end
end

subplot(311)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Tail beating')

subplot(312)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the swimming normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Swimming')

subplot(313)
xlim([-1.5 1.5])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Freq half cycle','FontSize',13)
title('Frequency')

sgtitle('MAX POINTS FOR EACH EPISODE')

saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
    strcat('Resume_max_point_episode')))
close all


 %%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY FOR ALL POINTS 
 % ---- threshold to all episodes inside the retrospective window
 clear
 user_settings

 threshold_calcium_response = 0.20;

 for iiii=1:length(list)

     nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

     [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
     CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
     sr=1/CALCIUMroiTS.deeplabcut.sr;


     % will be consider peaks of the calcium waves every wave that is at least
     % 20% amplitude of the bigger amplitud in the recording
     sum=0;  % for counting 1,2,3... in every loop independing of the index

     for i = CALCIUM.index_order % neurons order from medial to lateral. This variable is below

         temp1 = fieldnames(CALCIUM.AMPLITUDE_Y);
         myfield1 = char(temp1(i)); % extracting the name of each neuron
         sum = length(i) + sum;

         for h=1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
             temp2 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
             myfield2 = char(temp2(h)); % extracting the name of each episode
             index=CALCIUM.RESPONSE.(myfield1)(h); % This index will tell us if the cell is not responding
             % to the current episode and also at wich CALCIUM transient is
             % couple the response

             if any(index)         %isfield(CALCIUM.WINDOW_Y,(myfield1)) % Alternative for the if

                 % ------------- Next, ploting is fot the tail beating

                 if isfield(CALCIUM.WINDOW_Y.(myfield1),(myfield2)) % if the field exits, plot the values of the amplitudes
                     %of the swimming included in the window of that peak
                     a=find(CALCIUM.WINDOW_Y.(myfield1).(myfield2));
                     c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                     if  c>threshold_calcium_response*CALCIUM.MAXCALCIUM_RESPONSE
                         b=(c.*ones(1,length(CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                         subplot(311)

                         % Ploting all the amplitude of each episode

                         scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)))',...
                             CALCIUM.WINDOW_Y.(myfield1).(myfield2)(a)./ CALCIUM.NMAX_TAILBEATING,80,b','filled'), hold on
                     else
                     end


                 end
                 clear a b c

                 % Next is for the swimming amplitude
                 if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
                     if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
                         a=find(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2));
                         c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                         if  c>threshold_calcium_response*CALCIUM.MAXCALCIUM_RESPONSE
                             b=(c.*ones(1,length(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                             subplot(312)

                             % Ploting all the amplitude of each episode
                             scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)))',...
                                 CALCIUM.WINDOW_Y_SWIM.(myfield1).(myfield2)(a)./  CALCIUM.NMAX_SWIMMING,80,b','filled'), hold on



                             % Next is for the swimming freq
                             if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))

                                 if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))

                                     a=find(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2));
                                     if isempty(a)
                                     else
                                         b=(c.*ones(1,length(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)(a))));

                                         subplot(313)
                                         % Ploting all the amplitude of each episode
                                         scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)(a)))',...
                                             CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1).(myfield2)',80,b','filled'), hold on

                                     end
                                 else
                                 end
                             else
                             end

                         end
                     end
                 end
             else
             end
         end
     end
 end

subplot(311)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Tail beating')

subplot(312)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the swimming normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Swimming')

subplot(313)
xlim([-1.5 1.5])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Freq half cycle','FontSize',13)
title('Frequency')

sgtitle('ALL POINTS FOR EACH EPISODE') 

 saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
                strcat('Resume_Threshold_episode')))
 close all


 %%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
 % ----  all episodes inside the swimming episode (NO WINDOWS)
 
clear
user_settings
for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;


    % will be consider peaks of the calcium waves every wave that is at least
    % 20% amplitude of the bigger amplitud in the recording
    sum=0;  % for counting 1,2,3... in every loop independing of the index

    for i = CALCIUM.index_order % neurons order from medial to lateral. This variable is below

        temp1 = fieldnames(CALCIUM.AMPLITUDE_Y);
        myfield1 = char(temp1(i)); % extracting the name of each neuron
        sum = length(i) + sum;

        for h=1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
            temp2 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
            myfield2 = char(temp2(h)); % extracting the name of each episode
            index=CALCIUM.RESPONSE.(myfield1)(h); % This index will tell us if the cell is not responding
            % to the current episode and also at wich CALCIUM transient is
            % couple the response

            if any(index)         %isfield(CALCIUM.WINDOW_Y,(myfield1)) % Alternative for the if

                % ------------- Next, ploting is for the tail beating

                if isfield(CALCIUM.WINDOW_Y.(myfield1),(myfield2)) % if the field exits, plot the values of the amplitudes
                    %of the swimming included in the window of that peak
                    a=CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield2);
                    c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                    b=(c.*ones(1,length(a))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                    subplot(311)

                    % Ploting all the amplitude of each episode

                    scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(a))',...
                       a./ CALCIUM.NMAX_TAILBEATING,40,b','filled'), hold on


                end
                clear a b c

                % Next is for the swimming amplitude
                if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
                    if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
                        a=CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2);
                        c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
                        b=(c.*ones(1,length(a))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
                        subplot(312)

                        % Ploting all the amplitude of each episode
                        scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(a))',...
                            a./  CALCIUM.NMAX_SWIMMING,40,b','filled'), hold on


                        % Next is for the swimming freq
                        if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))

                            if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))

                                a=CALCIUM.FREQ_SWIM_Y.(myfield2);
                                if isempty(a)
                                else
                                    b=(c.*ones(1,length(a)));

                                    subplot(313)
                                    % Ploting all the amplitude of each episode
                                    scatter(CALCIUM.Normalized_x_distance(sum)*ones(1,length(a))',...
                                        a',40,b','filled'), hold on
                                end

                            else
                            end
                        else
                        end

                    end
                end
            else
            end
        end
    end
end

subplot(311)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Tail beating')

subplot(312)
xlim([-1.5 1.5])
ylim([0 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the swimming normalized','FontSize',13)
a=colorbar;
ylabel(a,'ΔF/F%','FontSize',16,'Rotation',270);
hColourbar.Label.Position(1) = 3;
title('Swimming')

subplot(313)
xlim([-1.5 1.5])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Freq half cycle','FontSize',13)
title('Frequency')

sgtitle('NO WINDOWS')

saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
    strcat('Resume_NO_WINDOW_episode')))
close all





%%    DISTRIBUTION (ACUMULATE FREQUENCY WITH BINS)

nbins = 7; % Select the number of bins that you want on your distribution 

color = ['b','g','b','k','c','m']; 
for i = 1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
    temp = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
    myfield = char(temp(i));
    temp3 (1,i) = min(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield));
    temp3 (2,i) = max(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield));
end 
temp1 = max(max(temp3)); 
temp2 = min(min(temp3)); 

vector = linspace(floor(temp2),ceil(temp1),nbins);


    for ii = 1:length(vector)-1
        creating_x_axis(ii) = (vector(ii) + vector(ii+1))/2;
    end


for i = 1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
    temp = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
    myfield = char(temp(i));






    creating_y_axis.(myfield) = NaN(length(vector)-1,...
        length(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)));

    for ii = 1:length(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield))

        for iii = 1:length(vector)-1
            if CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)(ii)>vector(iii) && ...
                    CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)(ii)<vector(iii+1)
                creating_y_axis.(myfield)(iii,ii) = 1;
            else
                creating_y_axis.(myfield)(iii,ii) = 0;

            end
        end


    end

  

    distribution = sum(creating_y_axis.(myfield),2);
    bar(creating_x_axis,distribution,color(i)), 
    alpha(0.2)
    hold on
%     plot(creating_x_axis,distribution,color(i),'LineWidth',4)

% area(creating_x_axis,distribution,'b')
% alpha(0.8)
xlim([0 70])



end

xlabel('Amplitude')
ylabel('#')













% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % FINAL PLOT
% 
% neuron_1 = 1;
% neuron_2 = 9;
% 
% figure
% subplot(611)
% plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x)
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylabel('x-coordinate')
% title('Tail Movement')
% subplot(612)
% %Normalize each Calcium wave using the max of each wave (optional, to have comparable sizes)
% % plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,neuron_1)./max(CALCIUMroiTS.diff_perc03.data(:,neuron_1 ))),hold on
% % plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,neuron_2)./max(CALCIUMroiTS.diff_perc03.data(:,neuron_2)))
% 
% plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,neuron_1)),hold on
% plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,neuron_2))
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylabel('%DelF/F')
% title('Calcium Imaging')
% subplot(613)
% scatter(middpoint(above_5),freq_sw(above_5),'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 0 0],...
%     'LineWidth',1.5)
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylim([0 40])
% ylabel('Hz')
% title('Instantaneus Frequency')
% 
% subplot(614)
% scatter(middpoint(above_5),real_amplitude,'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 0 0],...
%     'LineWidth',1.5)
% 
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylim([0 80])
% ylabel('x displacement')
% title('Instantaenus Amplitude')
% 
% 
% subplot(6, 1 ,[5 6])
% h = imagesc(perloc); colormap turbo; hold on
% xticks([])
% colorbar ('southoutside')
% 
% saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Final_figure_swimming',char(CALCIUM.list(nfish,1)))))
% 
% 
% 
% 
% 
% 
% 
% 
% % FINAL PLOT
% 
% figure
% subplot(611)
% plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x)
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylabel('x-coordinate')
% title('Tail Movement')
% subplot(612)
% %Normalize each Calcium wave using the max of each wave (optional, to have comparable sizes)
% % plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,1)./max(CALCIUMroiTS.diff_perc03.data(:,1))),hold on
% % plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,9)./max(CALCIUMroiTS.diff_perc03.data(:,9)))
% 
% plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,neuron_1)),hold on
% plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,neuron_2))
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylabel('%DelF/F')
% title('Calcium Imaging')
% subplot(613)
% scatter(middpoint,freq_sw,'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 0 0],...
%     'LineWidth',1.5)
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylim([0 40])
% ylabel('Hz')
% title('Instantaneus Frequency')
% 
% subplot(614)
% scatter(middpoint,real_amplitude,'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceColor',[1 0 0],...
%     'LineWidth',1.5)
% 
% xlim([CALCIUMroiTS.diff_perc03.times(1) CALCIUMroiTS.diff_perc03.times(end)])
% ylabel('x displacement')
% title('Instantaenus Amplitude')
% 
% 
% subplot(6, 1 ,[5 6])
% h = imagesc(perloc); colormap bone; hold on
% xticks([])
% colorbar ('southoutside')
% 
% saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Final_figure_tail_beating',char(CALCIUM.list(nfish,1)))))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%% RASTER 
clear
user_settings
figure('units','normalized','outerposition',[0 0 1 1])
for iiii=1:length(list)

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;

    imagesc(perloc./max(perloc,[],2))
for i=1:length(fieldnames(CALCIUM.EPISODE_X))
    temp = fieldnames(CALCIUM.EPISODE_X);
    myfield = char(temp(i)); 
    temp_xline = (CALCIUM.EPISODE_X.(myfield)(1)*1/CALCIUMroiTS.deeplabcut.sr)/CALCIUM.srate;
    xline(temp_xline,'k','LineWidth',3)
     temp_xline = (CALCIUM.EPISODE_X.(myfield)(end)*1/CALCIUMroiTS.deeplabcut.sr)/CALCIUM.srate;
    xline(temp_xline,'r','LineWidth',3)

end 

yline(CALCIUM.midline_respect_neurons,'white','LineWidth',3);
set(gca,'YTick',1:length(CALCIUM.index_order));% something like this
set(gca,'YTickLabels',[CALCIUM.index_order]);% something like this
colormap('parula')

saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Raster_cells',char(CALCIUM.list(nfish,1)),'.fig')))



%% PLOTING PHASIC AND ANTIPHASIC ACTIVITY
figure('units','normalized','outerposition',[0 0 1 1])
% Left hemisphere

%     Antiphasic
antiphasic_left = []; 

subplot(341)
plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,antiphasic_left),'color',[0.8500 0.3250 0.0980]),hold on
ylabel('% \Delta F')
yyaxis right
plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x,'color',[0 0.4470 0.7410])
title ('Antiphasic LH neurons')
ylabel('\Deltax Coordinates tail3 [px]')

%     Phasic

phasic_left = [15 16 19 20 24 25 26 30 34 ] ; 

subplot(345)
plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,phasic_left),'color',[0.4660 0.6740 0.1880]),hold on
ylabel('% \Delta F')
yyaxis right
plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x,'color',[0 0.4470 0.7410])
title ('Phasic LH neurons')
ylabel('\Deltax Coordinates tail3 [px]')


subplot(349)
plot(CALCIUMroiTS.diff_perc03.times,mean(-CALCIUMroiTS.diff_perc03.data(:,antiphasic_left),2),'color',[0.8500 0.3250 0.0980]),hold on
ylabel('% \Delta F')
plot(CALCIUMroiTS.diff_perc03.times,mean(-CALCIUMroiTS.diff_perc03.data(:,phasic_left),2),'color',[0.4660 0.6740 0.1880])
yyaxis right
plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x,'color',[0 0.4470 0.7410])
ylabel('\Deltax Coordinates tail3 [px]')
title ('Average Phasic & Anphasic LH neurons')






%% Right hemisphere

%     Antiphasic
antiphasic_right = [1]; 

subplot(342)
plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,antiphasic_right),'color',[0.8500 0.3250 0.0980]),hold on
ylabel('% \Delta F')
yyaxis right
plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x,'color',[0 0.4470 0.7410])
title ('Antiphasic RH neurons')
ylabel('\Deltax Coordinates tail3 [px]')


% Right hemisphere

phasic_right = [2 3 4 5 6]; 


subplot(346)
plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(:,phasic_right),'color',[0.4660 0.6740 0.1880]),hold on
ylabel('% \Delta F')
yyaxis right
plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x,'color',[0 0.4470 0.7410])
title ('Phasic RH neurons')
ylabel('\Deltax Coordinates tail3 [px]')


subplot(3,4,10)
plot(CALCIUMroiTS.diff_perc03.times,mean(-CALCIUMroiTS.diff_perc03.data(:,antiphasic_right),2),'color',[0.8500 0.3250 0.0980]),hold on
ylabel('% \Delta F')
plot(CALCIUMroiTS.diff_perc03.times,mean(-CALCIUMroiTS.diff_perc03.data(:,phasic_right),2),'color',[0.4660 0.6740 0.1880])
yyaxis right
plot(CALCIUMroiTS.deeplabcut.lasertimesync,CALCIUMroiTS.deeplabcut.tail3.x,'color',[0 0.4470 0.7410])
ylabel('\Deltax Coordinates tail3 [px]')
title ('Average Phasic & Anphasic RH neurons')


% Ploting the position of antiphasic


subplot(3,4,[3 4 7 8 11 12])
% Antiphasic 
temp_vector = [antiphasic_left'; antiphasic_right']; 
imagesc(CALCIUM.mean), hold on
for i = temp_vector'
temp = CALCIUM.roi.manual_poly{i}; 
c = [0.8500 0.3250 0.0980]; 
fill(temp(:,1),temp(:,2),c),
end 

% Ploting the position of antiphasic


% phasic 
temp_vector = [phasic_left'; phasic_right']; 
for i = temp_vector'
temp = CALCIUM.roi.manual_poly{i}; 
c = [0.4660 0.6740 0.1880]; 
fill(temp(:,1),temp(:,2),c),
end 
sgtitle(char(CALCIUM.list(nfish,1)))
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Phasic_Antiphasic_',char(CALCIUM.list(nfish,1)),'.fig')))
close all

CALCIUM.Phase.phasic_left = phasic_left; 
CALCIUM.Phase.phasic_right = phasic_right;
CALCIUM.Phase.antiphasic_left = antiphasic_left; 
CALCIUM.Phase.antiphasic_right = antiphasic_right; 
CALCIUMimg('save',CALCIUM,[],list,nfish);


