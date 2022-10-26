function [] = roi_preview(backgr, coord)

% coord: coordinates matrix (not cell structure)

figure
imagesc(backgr); colormap('bone'); hold on

    plot(coord(:,1), coord(:,2), 'color',[.1 .5 .5], 'LineWidth', 1.5); hold on

end

%% Created: 06/02/21 (from polygon_preview)
% Updated 