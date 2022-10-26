function [roi_idx] =roi_idxfromname(VSDI,roi_names)

% INPUT
% 'roi_names' - 
% OUTPUT
% 'roi_idx'

for roii= 1:length(roi_names) 
roi_idx = strcampi(roi_names{roii},) ;
end
end