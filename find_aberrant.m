function [reject_idx]  = find_aberrant(GS, cond_idx, plot_on)

% Find idx of trials to reject based on the criterion of Global Signal
% Substraction.

% INPUT: 
%  'GS': 2d array (timepoints, trials). Contains the whole-brain
%   pixels-average (GS) for each trial
% 'cond_idx': indexes of the condition. The GS substraction is made
% condition-wise.
% 'plot_on': 1 if we want to plot the result

% substract from 

% STEP 0 Preprocess and crop brain (p5)
% The function: 
    % STEP 1  Average GS across trials from one condition 
    % STEP 2 (for each trial): substract the averaged GS and sum the
    % residuals
% % % % %     
% % % % %     for trialidx = trialset
% % % % % cropmovie = medianmovie.*cropmask3d;
% % % % % movies_p8(:,:,:,trialidx) = cropmovie;
% % % % % 
% % % % % % compute GS for that movie at each timepoint
% % % % % for ti = 1:t
% % % % % cropframe = movie4d(:,:,ti).*cropmask;
% % % % % % GS(ti,trialidx) =  sum(cropframe(:))/sum(cropmask(:));
% % % % % GSbase(ti,
% % % % % 
% % % % % end
% % % % % 
% % % % % clear roimask rasterclus_p8
% % % % %   end %end loop through roi
    
    % STEP 3 : compute the mean and Std of the accumulated residuals 
    % STEP 4: Discard trials if the accumulated residuals for that trial is superior to
    % mean + 2*Std 

% input control
if nargin <3
    plot_on = 0;
end

% If condition index is a column, turn into row to avoid problems with the
% loop
if size(cond_idx,1) >  size(cond_idx,2) 
cond_idx = cond_idx'; 
end
% end input control

% get GS from selected trials
GSave = mean(GS(:,cond_idx),2); %plot(GSave) 
std_crit = 2; 

counti = 1; 
for triali = cond_idx
resid(:,triali) = GS(:,triali) - GSave;
accum_resid(triali) = sum(abs(resid(:,triali))); 
end

% to calculate the threshold, just idx of the trials should be included,
% otherwise, 0-values from other trials will be computed wrongly into the
% mean and std
thresh = mean(accum_resid(cond_idx)) + std_crit*std(accum_resid(cond_idx)); %threshold for rejection

reject_idx =[];
for triali= cond_idx
   if accum_resid(triali) > thresh
       reject_idx = [reject_idx triali]; 
   else
       reject_idx = reject_idx; 
   end
end

if plot_on
    figure
    % plot rejected vs non-rejected
    cleanidx = setdiff(cond_idx, reject_idx); 
    plot(GS(:,cleanidx), 'color', '#F3E5AB'); hold on;
    plot(resid(:,reject_idx), 'color', '#97978F'); hold off
    % title(strcat('n reject=', num2str(length(reject_idx))))
end

end

%% Updated: 26/10/20