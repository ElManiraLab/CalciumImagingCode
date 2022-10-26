function [output] = SingleQuotation(input)
%SINGLEQUOTATION   [SINGLE_QUOTATION] = SINGLE_QUOTATION [DOUBLE_QUOTATION]
%   Elimina las comillas dobles y mantiene las simples. Avisa de algún
%   error en la ruta
if  size(find(input=='"')) == [1 0] % Implica que no existe doble comilla (Vector vacio)
    output=input; % Por tanto existe comilla simple y se mantiene la notación
else   % Implica que sí existe doble comilla
    Posicion = find(input=='"'); % Genera un output de dónde están las dobles comillas
    if length(Posicion)==1 %Si solo hay una comilla
        if  Posicion(1)==1 % Está al principio
            output = input(Posicion(1)+1:length(input));
        elseif Posicion(1)==length(input) % Está al final
            output = input(1:length(input)-1);
        else
            error('Revisar ruta: doble comilla dentro de la ruta') % Este entre medias
        end
        
    elseif length(Posicion)==2
        if  Posicion(1)~=1 | Posicion(2)~=length(input) %Solo si no están al principio o final lanza error
            error('Revisar ruta')
        elseif  Posicion(1)==1 & Posicion(2)==length(input) %Solo si están al principio o final ejecuta la linea
            output = input(Posicion(1)+1:Posicion(2)-1);
        else
            error('Revisar ruta') % Si nada de lo anterior es cierto lanza errror
        end
    else 
        error('Revisar ruta: más de dos dobles comillas') % Si hay más de dos dobles comillas
    end
end
end


