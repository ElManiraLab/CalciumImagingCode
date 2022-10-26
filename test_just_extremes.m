%%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
 % ----  max pf each episodes 
 % ----  (NO WINDOWS)
 
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

    for i = [1 2] % neurons order from medial to lateral. This variable is below

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

                    scatter(CALCIUM.Normalized_x_distance(sum),...
                       max(a)./ CALCIUM.NMAX_TAILBEATING,40,b(1),'filled'), hold on


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
                        scatter(CALCIUM.Normalized_x_distance(sum),...
                            max(a)./  CALCIUM.NMAX_SWIMMING,40,b(1),'filled'), hold on


                        % Next is for the swimming freq
                        if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))

                            if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))

                                a=CALCIUM.FREQ_SWIM_Y.(myfield2);
                                if isempty(a)
                                else
                                    b=(c.*ones(1,length(a)));

                                    subplot(313)
                                    % Ploting all the amplitude of each episode
                                    scatter(CALCIUM.Normalized_x_distance(sum),...
                                        max(a),40,b(1),'filled'), hold on
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


%%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
 % ----  max pf each episodes 
 % ----  (NO WINDOWS)
 % --- normalized to each cell the calcium transient
 
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

                    scatter(CALCIUM.Normalized_x_distance(sum),...
                       max(a)./ CALCIUM.NMAX_TAILBEATING,150,b(1)./max(CALCIUM.BASELINE_TO_PEAK.(myfield1)),'filled'), hold on


                end
                clear a b c

                % Next is for the swimming amplitude
                
%                 if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
%                     if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
%                         a=CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2);
%                         c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
%                         b=(c.*ones(1,length(a))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
%                         subplot(312)
% 
%                         % Ploting all the amplitude of each episode
%                         scatter(CALCIUM.Normalized_x_distance(sum),...
%                             max(a)./  CALCIUM.NMAX_SWIMMING,150,b(1),'filled'), hold on
% 
% 
%                         % Next is for the swimming freq
%                         if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))
% 
%                             if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))
% 
%                                 a=CALCIUM.FREQ_SWIM_Y.(myfield2);
%                                 if isempty(a)
%                                 else
%                                     b=(c.*ones(1,length(a)));
% 
%                                     subplot(313)
%                                     % Ploting all the amplitude of each episode
%                                     scatter(CALCIUM.Normalized_x_distance(sum),...
%                                         max(a),150,b(1),'filled'), hold on
%                                 end
% 
%                             else
%                             end
%                         else
%                         end
% 
%                     end
%                 end
            else
            end
        end
    end
end

subplot(311)
xlim([-1.5 1.5])
ylim([0.7 1])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
colormap('turbo')
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


%%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
 % ----  max pf each episodes 
 % ----  (NO WINDOWS)
 % --- normalized to each cell the calcium transient
 % ---- plotting the no response like empty circle

 
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
                   

                    % Ploting all the amplitude of each episode

                    scatter(CALCIUM.Normalized_x_distance(sum),...
                       max(a)./ CALCIUM.NMAX_TAILBEATING,150,b(1)./max(CALCIUM.BASELINE_TO_PEAK.(myfield1)),'filled'), hold on


                end
                clear a b c

                % Next is for the swimming amplitude
%                 if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
%                     if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
%                         a=CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2);
%                         c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
%                         b=(c.*ones(1,length(a))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
%                         subplot(312)
% 
%                         % Ploting all the amplitude of each episode
%                         scatter(CALCIUM.Normalized_x_distance(sum),...
%                             max(a)./  CALCIUM.NMAX_SWIMMING,150,b(1),'filled'), hold on
% 
% 
%                         % Next is for the swimming freq
%                         if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))
% 
%                             if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))
% 
%                                 a=CALCIUM.FREQ_SWIM_Y.(myfield2);
%                                 if isempty(a)
%                                 else
%                                     b=(c.*ones(1,length(a)));
% 
%                                     subplot(313)
%                                     % Ploting all the amplitude of each episode
%                                     scatter(CALCIUM.Normalized_x_distance(sum),...
%                                         max(a),150,b(1),'filled'), hold on
%                                 end
% 
%                             else
%                             end
%                         else
%                         end
% 
%                     end
%                 end
            else
                 
                a=CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2);
                scatter(CALCIUM.Normalized_x_distance(sum),...
                       max(a)./ CALCIUM.NMAX_TAILBEATING,150,[0.3010 0.7450 0.9330]), hold on

            end
        end
    end
end

subplot(311)
xlim([-1.5 1.5])
ylim([0.8 1.02])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
colormap('turbo')
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



%%  PLOTING AMPLITUDE TAIL BEATING, AMPLITUD SWIMMING AND FREQUENCY 
 % ----  max pf each episodes 
 % ----  (NO WINDOWS)
 % --- normalized to each cell the calcium transient
 % ---- plotting the no response like empty circle
  % ---- just one side
 
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
                   

                    % Ploting all the amplitude of each episode

                    scatter(abs(CALCIUM.Normalized_x_distance(sum)),...
                       max(a)./ 58.48,150,b(1)./max(CALCIUM.BASELINE_TO_PEAK.(myfield1)),'r','filled'), hold on


                end
                clear a b c

                % Next is for the swimming amplitude
%                 if isfield(CALCIUM.WINDOW_Y_SWIM,(myfield1))
%                     if isfield(CALCIUM.WINDOW_Y_SWIM.(myfield1),(myfield2))
%                         a=CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2);
%                         c=CALCIUM.BASELINE_TO_PEAK.(myfield1)(index);
%                         b=(c.*ones(1,length(a))); %./CALCIUM.MAXCALCIUM_RESPONSE_TOTAL; %Normalize
%                         subplot(312)
% 
%                         % Ploting all the amplitude of each episode
%                         scatter(CALCIUM.Normalized_x_distance(sum),...
%                             max(a)./  CALCIUM.NMAX_SWIMMING,150,b(1),'filled'), hold on
% 
% 
%                         % Next is for the swimming freq
%                         if isfield(CALCIUM.WINDOW_X_SWIM_FREQ,(myfield1))
% 
%                             if isfield(CALCIUM.WINDOW_X_SWIM_FREQ.(myfield1),(myfield2))
% 
%                                 a=CALCIUM.FREQ_SWIM_Y.(myfield2);
%                                 if isempty(a)
%                                 else
%                                     b=(c.*ones(1,length(a)));
% 
%                                     subplot(313)
%                                     % Ploting all the amplitude of each episode
%                                     scatter(CALCIUM.Normalized_x_distance(sum),...
%                                         max(a),150,b(1),'filled'), hold on
%                                 end
% 
%                             else
%                             end
%                         else
%                         end
% 
%                     end
%                 end
            else
                 
                a=CALCIUM.FINAL_POINTS_SWIMMING_Y.(myfield2);
                scatter(abs(CALCIUM.Normalized_x_distance(sum)),...
                       max(a)./ 58.48,150,'b','filled'), hold on

            end
        end
    end
end

subplot(311)
xlim([-1.5 1.5])
ylim([0.8 1.02])
xlabel('medial (0) to lateral(n) (# neurons)','FontSize',13)
ylabel('Amplitude of the tail beating normalized','FontSize',13)
colormap('turbo')
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