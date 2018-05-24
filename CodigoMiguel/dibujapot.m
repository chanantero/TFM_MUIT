function [ p_final ] = dibujapot( SIGTOTAL,L, alt, ad_sources, mallado_x, mallado_y,vent,desp,fileavi);
%dibujapot
% SIGTOTAL -> matriz con los datos de la señal en cada columna
% L -> Tamaño de la sala [x,y,z];
%   alt -> coordenadas de los altavoces
%   ad_sources -> coordenadas de las fuentes adicionales
%   mallado_x -> vector con las posiciones de los puntos de control en la coordenada X 
%   mallado_y -> vector con las posiciones de los puntos de control en la
%       coordenada Y
%   Detailed explanation goes here

[num_muestras num_mic]=size(SIGTOTAL);

p=zeros(length(mallado_x),length(mallado_y));
PSIGTOTAL=SIGTOTAL.^2;


maxi=10*log10(max(max(PSIGTOTAL)));
mini=10*log10(min(min(PSIGTOTAL(num_muestras-vent:num_muestras,:))));
mini=max([mini maxi-50]);

scrsz = get(0,'ScreenSize');
fig=figure('Position',[10 61 scrsz(3)-20 scrsz(4)-180]);
plot(alt(1,1:96),alt(2,1:96),'*');
axis([0 L(1) 0 L(2)]);
if exist('ad_sources')==0
    else 
    hold on
    %Altavoces adicionales 
    plot(ad_sources(:,1),ad_sources(:,2),'o');
end

if (exist('fileavi')==1);
%fig=figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
set(fig,'DoubleBuffer','on');
mov = avifile([fileavi '.avi']);
end


num_pot=1+fix((num_muestras-vent)/desp);

for i=1:num_pot
pot=sum(PSIGTOTAL(1+desp*(i-1):vent+desp*(i-1),:))./vent;
%pot(:,1:length(mallado_x))=1;
p=reshape(10*log10(pot),length(mallado_x),length(mallado_y))';

pcolor(mallado_x,mallado_y,p); shading interp; lighting phong; colormap hot; axis off; %zlim([0 1]); 
caxis([mini maxi]);
colorbar;
%set(gcf,'Color', [0 0 0], 'Number', 'off', 'Name', sprintf('FDTD, iteracion = %i', n));
title(sprintf('Muestra = %.2f',vent+(i-1)*desp),'Color',[1 0 0],'FontSize', 22); 
drawnow;

if (exist('fileavi')==1);
F = getframe(gca);
       	mov = addframe(mov,F);
end
end

pot=sum(PSIGTOTAL(end-vent:end,:))./vent;
p=reshape(10*log10(pot),length(mallado_x),length(mallado_y))';
pcolor(mallado_x,mallado_y,p); shading interp; lighting phong; colormap hot; axis off; zlim([0 1]); 
caxis([mini maxi]);
colorbar;
%set(gcf,'Color', [0 0 0], 'Number', 'off', 'Name', sprintf('FDTD, iteracion = %i', n));
title(sprintf('Muestra = %.2f',num_muestras),'Color',[1 0 0],'FontSize', 22); 
drawnow;

if (exist('fileavi')==1);
   mov = close(mov); 
end

p_final=p;
end

