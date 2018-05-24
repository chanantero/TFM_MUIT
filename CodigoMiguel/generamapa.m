function [POT_ad,POT_ar]=generamapa(xsignal,h_array,filtros_array,activo_array,xrefsignal,h_ad_sources);

% Genera los valores de la potencia instantanea en todos los puntos de
% control de las se�ale definidas en xsignal y xrefsignals  
% La salida ser� una matriz con tantas filas como duraci�n de la se�al y
% tantas columnas como puntos de control
%
% xsignal-> vector columna con del tama�o de la se�al (ser� sintetizado por el array)
% h_array .> estructura de datos con los caminos ac�sticos de los altavoces
%   del array hasta los puntos de control
%   filtros_Array -> Matriz con los filtros del array
%   vector para indicar los elementos activos del array
% xrefsignal-> vector con tantas columnas como funtes adicionales y tantas filas como el tama�o de la se�al
% h_ad_sources -> estructura de datos con los caminos ac�sticos de las
% fuentes adicionales hasta los puntos de control
% 





[durref,num_sec_ad]=size(xrefsignal);
dur=min([durref length(xsignal)]);

[L,num_alt_ad,nx1,ny1]=size(h_ad_sources);
[L,num_alt_arr,nx2,ny2]=size(h_array);

nx=max([nx1 nx2]);
ny=max([ny1 ny2]);


POT_ad=zeros(dur,nx*ny);
POT_ar=zeros(dur,nx*ny);

for altavoz=1:num_alt_ad
    cont=1;
    for y=1:ny
        for x=1:nx
        POT_ad(:,cont)=POT_ad(:,cont)+filter(h_ad_sources(:,altavoz,x,y),1,xrefsignal(:,altavoz));
        cont=cont+1;
        end
    end
end


for altavoz=1:length(activo_array)
    cont=1;
    xfsignal=filter(filtros_array(:,activo_array(altavoz)),1,xsignal);
    for y=1:ny
        for x=1:nx
        POT_ar(:,cont)=POT_ar(:,cont)+filter(h_array(:,activo_array(altavoz),x,y),1,xfsignal);
        cont=cont+1;
        end
    end
end



