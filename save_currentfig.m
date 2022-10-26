function save_currentfig (savepath, name)
% save_cfig (savepath, name) saves the opened figures in 'savepath' and
% closes them

h =  findobj('type','figure');
n = length(h);

for ii = 1:n
   name = strcat(name, num2str(h(ii).Number)); 
   saveas(h(ii), [savepath, name, 'jpg'],'jpg')
   close(h(ii))
end

%% Created: 07/02/21
