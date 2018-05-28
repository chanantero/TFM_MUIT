
fs = 8000;
c=340;

fte = [0.1 2.3 1.65;0.1 1.3 1.65;0.1 3.3 1.65];
%fte=[alt(1,93)-0.001 alt(2,93) 1.65;2,2.3,1.65];
%fte=alt(:,93);
fuente_ruido=1;

% Definimos la sala y las fuentes de ruido
disp('Calculando las respuestas del sistema acústico a modelar...')
% [h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,0.05,0.05,0,250,fs,c,fte,1);
% save salasinrev h_array h_ad_sources alt tecta fs fte L mallado_x mallado_y;
% 
% 
% [h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,0.2,0.2,0.15,250,fs,c,fte);
% save salaconrev2 h_array h_ad_sources alt tecta fs fte L mallado_x mallado_y;
[h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,0.2,0.2,0,250,fs,c,fte);
%save salasinr3 h_array h_ad_sources alt tecta fs fte L mallado_x mallado_y;
% Calculamos las funciones directoras del array definido por 'alt' para sintentizar una fuente en la posición fte. 


disp('Calculando la configuración del array...')
[filtros_array,an,tn,activo_array]=WFS_DrivingSignals(alt,fte(fuente_ruido,:),c,fs);

% Definimoa la señal de ruido
x=sin(2*pi*0.01*[0:1999])';
%x=randn(2000,1);
realFreq = 0.01*fs;

% Definimos la señal de la fuente virtual;
xv=-0.0957*x;
%xv=-x;
xv(1:999)=0; % empieza a actuar en la muestra 1000;
%x=[zeros(30,1);x(1:end-30)];



% Calculamos al potencia en los puntos de control
disp('Calculando las señales en los puntos de control...')
[POT_ad,POT_ar]=generamapa(xv,h_array,filtros_array,activo_array,x,h_ad_sources(:,fuente_ruido,:,:));

% Representamos el mapa de presiones
 dv=50; % Duración de la ventana para promediar al potencia
 despv=25;  % Desplazamiento del enventanado para realizar el calculo de la
% potencia

disp('Dibujando el mapa de potencias');
dibujapot(POT_ad+POT_ar ,L, alt, fte(fuente_ruido,:), mallado_x, mallado_y,dv,despv);

% 
numSamp = size(POT_ad, 1);
t = (0:numSamp - 1)/fs;
ax = axes(figure);
plot(ax, t, POT_ad(:, 1), t, POT_ar(:, 1))