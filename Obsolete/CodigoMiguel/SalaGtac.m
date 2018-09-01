function [h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources]=SalaGtac(z,delta_x,delta_y,t_r,Np,fs,c,ad_sources,ok);
%[h_array,h_ad_sources ] = SalaGtac(z,delta_x,delta_y,t_r,Np,fs,ad_sources);
%  Función que devuelve las respuestas al impulso de los altavoces de la sala
%  GTAC (y opcionalmente otros altavoces adicionales) respecto de un mallado de 
%  puntos de control definidos en un plano a la altura z y equiespaciados delta_x
% y delta_y 
%   También devuelve las coordenadas y la orientación de los altavoces en
%       las matrices alt y tecta, así como las dimensiones de la sala en la
%       variable L
%       En las variables mallado_x y mallado_y se devuelven los puntos
%       donde se quiere monitorizar la señal (hacia donde se miden las
%       respuestas al impulso).
%   Las respuestas al impulso se almacenan en la forma
%   h_array(:,altavoz,posx,posy)
%
%       z-> Altura del mallado de control (z=1.65 para colocarlo a la altura de
%           los altavoces).
%       delta_x -> distancia entre puntos de control en la dirección x (0.2 para simular las medidas de Kristoff)
%       delta_y -> distancia entre puntos de control en la dirección y (0.2 para simular las medidas de Kristoff)
%       t_r -> tiempo de reverberación de la sala (0 para simular espacio
%           libre).
%       Np -> Nº de puntos para las respuestas al impulso 
%       fs -> Frecuencia de trabajo (por defecto 44100)
%       c -> velocidad del sonido (por defecto 340)
%       ad_sources -> Matriz de tamaño Kx3 con las coordenadas x,y,z de las posiciones de
%           las K fuentes adicionales. 
%  Si no se desea calcular las respuestas al impulso del array sino sólo de
%  las fuentes adicionales, la variable ok debe ser no nula


if exist('fs')==0
    fs=44100;
    disp('fs (frecuencia de muestreo)=44100')
end

if exist('c')==0
    c=340;
    disp('c (velocidad del sonido)=340 m/s')
end

z_sala=2.64; %Altura de la sala
x_sala=9.13; %Dimesión larga de la sala
y_sala=4.48; %Dimesión corta de la sala


posAlt1=[x_sala-0.56-6.36 y_sala-0.46-1.0182]; % Posición del altavoz de referencia 1




alt_z=1.65; % Coordenada z (altura) de los altavoces del array
% Coordenadas x,y de los altavoces del array
alt(2,1:8)=(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,1:8)=(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,9:32)=1.0182;
alt(1,9:32)=(1.0182)+(0.09:0.18:2*1.44+1.35);
alt(2,33:40)=-(0.09:0.18:1.35)*cos(45*pi/180)+1.0182;
alt(1,33:40)=(1.0182+3*1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,41:48)=-(0.09:0.18:1.35);
alt(1,41:48)=(2*1.0182+3*1.44);
alt(2,49:56)=-(1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,49:56)=(2*1.0182+3*1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,57:80)=-(1.0182+1.44);
alt(1,57:80)=(1.0182)+(2*1.44+1.35:-0.18:0.09);
alt(2,81:88)=-(1.0182+1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,81:88)=1.0182-(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,89:96)=-(1.44)+(0.09:0.18:1.35);
alt(1,89:96)=0;

alt(1,:)=alt(1,:)+posAlt1(1);
alt(2,:)=alt(2,:)+posAlt1(2);
% Orientación de los altavoces
tecta=[ones(1,8)*45,ones(1,24)*0,-ones(1,8)*45,-ones(1,8)*90,ones(1,8)*-135,ones(1,24)*-180,-ones(1,8)*225,ones(1,8)*90];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Definición del mallado
Lim_inf_izq=alt(:,82)+[0.1;0.1];
Lim_sup_der=alt(:,34)-[0.1;0.1];


mallado_x=Lim_inf_izq(1):delta_x:Lim_sup_der(1);
mallado_y=Lim_inf_izq(2):delta_y:Lim_sup_der(2);


% Lim_inf_izq=alt(:,82)
% Lim_sup_der=alt(:,34)
% mallado_x=(Lim_inf_izq(1)+Lim_sup_der(1))/2;
% mallado_y=(Lim_inf_izq(2)+Lim_sup_der(2))/2;



nx=length(mallado_x);
ny=length(mallado_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Representación gráfica

%Altavoces array
figure,plot(alt(1,1:96),alt(2,1:96),'*');
axis([0 x_sala 0 y_sala]);
if exist('ad_sources')==0
    disp('No hay altavoces adicionales')
else 
    hold on
    %Altavoces adicionales 
    plot(ad_sources(:,1),ad_sources(:,2),'o');
end


hold on
% Puntos de control
 Mnx=repmat(mallado_x,1,ny);
 Mny=repmat(mallado_y,nx,1);
 Mny=Mny(:)';
 plot(Mnx,Mny,'.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cálculo de las respuestas al impulso
%c=340;
L=[x_sala,y_sala,z_sala];
beta=t_r;

r=zeros(nx*ny,3);
pos=0;

for y=1:ny
    for x=1:nx
    pos=pos+1;
    r(pos,:)=[mallado_x(x),mallado_y(y),z];
    end
end

% para las del array:

if exist('ok')==1
    disp('Sólo se calculan las respuestas al impulso de los altavoces adicionales')
h_array=0;
else
h_array=zeros(Np,96,nx,ny);

 for a=1:96
 s=[alt(1,a),alt(2,a),alt_z]; 
 h = rir_generator(c, fs, r, s, L, beta, Np)';
 h_array(:,a,:,:)=reshape(h,Np,nx,ny);
 end
end


if exist('ad_sources')==0
   % disp('No hay altavoces adicionales')
else 
    %Para los altavoces adicionales 
    
    [num_alt_ad,dim]=size(ad_sources);
    h_ad_sources=zeros(Np,num_alt_ad,nx,ny);
    for a=1:num_alt_ad
        s=ad_sources(a,:);
    h = rir_generator(c, fs, r, s, L, beta, Np)';
    h_ad_sources(:,a,:,:)=reshape(h,Np,nx,ny);
    end
end




end

