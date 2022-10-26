function rowvector = makeRow(vector)
% Checks whether vector is a column and, if so, turns it into a row

% input control
if size(size(vector))> 2
    error('input a 2d vector')
end
% end

if length(vector) >1 % If it is a scalar, it leaves it unchanged

   if size(vector,1) > size(vector,2)
      rowvector = vector';
   else 
       rowvector = vector;
   end
   
else
       rowvector = vector;

end

end