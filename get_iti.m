function [iti] = get_iti(VSDI)
% [iti] = get_iti(VSDI) gets, for each trial, the previous ITI in seconds

% INPUT: 'VSDI' structure with field 'VSDI.list'; whose field 'Date'
% contains the datetime for each trial

ntri = length(VSDI.list);
for t= 1:ntri
   tiempos(t) = VSDI.list(t).Date; 
end
tiempos = tiempos';

difsec = seconds(diff(tiempos));
iti = NaN(size(tiempos)); 
iti(2:end) = difsec;
iti=fillmissing(iti, 'constant',Inf);%insert first value
% convertir a segundos

end

