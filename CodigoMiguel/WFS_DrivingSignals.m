function [pout_sign_aux,an,tn,parray_act]=WFS_DrivingSignals2(alt,fte,c,fs);


%%%% El objetivo de este código es calcular las amplitudes y los retardos
%%%% de cada altavoz de un array en configuración octogonal de 96 altavoces
%%%% cuyas posiciones están definidas en la matriz alt de dos filas y 96
%%%% columnas (una fila para cada coordenada cartesiana y una columna para
%%%% cada altavoz) dada una posición de una fuente virtual en la variable
%%%% fte (que contiene las coordenadas x e y de la posición de la fuente). 
%%%% 
%%%% Además, el código calcula por trigonomtría qué altavoces deberían
%%%% estar activos en función del angulo relativo de incidencia del sonido
%%%% procedeente de la fuente virtual. 
%%%%
%%%% Las variables de salida son:
%%%% parray_act -> contiene el identificador de los altavoces activos. Ej:
%%%% si todos estuviesen activos parra_act=[1:96];
%%%%
%%%% an-> un vector con las amplitudes con las que habría que ponderar la
%%%% señal de la fuente virtual en cada altavoz.
%%%% tn-> un vector con los retardos que habría que aplicar a la señal de
%%%% la fuente virtual en cada altavoz. 
%%%%
%%%% Por último, en pout_sign_aux (matriz de 512 filas y tantas columnas
%%%% como altavoces tenga el array) se almacena un filtro fir de 512
%%%% elementos para cada altavoz con la información de la amplitud y
%%%% retardo que hemos calculado anteriormente. De esta forma podemos usar
%%%% estos filtros para filtrar la señal de la fuente virtual antes de
%%%% emitirla por cada altavoz. 

%%%% El código, por sí mismo, no maneja señales. Es decir, no crea la
%%%% sintesis del campo sonoro. Sólo calcula los factores de corrección (amplitud y retardo)que
%%%% habría que aplicarle a la señal para generala por cada altavoz y
%%%% obtener la sintesis de la fuente virtual. Estos factores sevirían
%%%% para, dada una señal, calcular las señales que debería generar cada
%%%% altavoz, pero esa parte no está implementada. Si se usan los filtros diseñados, simplemente habría que filtrar la señal por cada filtro antes de emitirla por el altavoz correspondiente.  



%c=343; % Velocidad del sonido
comp=sqrt(-1); % j
ang_max=90*(pi/180); % Angulo máximo de incidencia del sonido (90º) en radianes
ord=60;
Ax=0.18; % Separación entre altavoces
%FS=44100; % Frecuencia de trabajo
land=-1; % porque está fuera del array de altavoces
ali_frec=c/(2*max(Ax)*sin(ang_max)); % frecuencia de alias
ali_frec=ceil(ali_frec); % frecuencia de alias redondeada al alza

% Diseño de filtro para la reconstrucción de la señal en campo cercano
% WFS podría funcionar sin este filtro
frec_filtro=linspace(0,ceil(fs/2)+1,(ceil(fs/2)/40)+1);
longit=length(frec_filtro)-mod(length(frec_filtro),2);
amp_filtro=sqrt(-(frec_filtro(1:ceil((ali_frec+1000)/40)))/(c*comp));
vector(1:ceil((ali_frec+1000)/40))=amp_filtro;
vector(ceil((ali_frec+1000)/40)+1:longit)=amp_filtro(end);
H_1=firls(ord,frec_filtro(1:longit)/frec_filtro(longit),abs(vector)); % H_1 contiene la respuesta al impulso del filtro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ret_H_1=ord/2;
%figure, plot(H_1)
%figure,freqz(H_1)

%%%%% Definición de la geometría del array (que se corresponde con el de la sala GTAC)
%% Coordenadas x,y de las posiciones de los altavoces (el eje x sería paralelo a la pared larga) 
%% x=0 se corresponde con la línea larga de los altavoces en la pared de enfrente de la puerta
%% y=0 se corresponde con la línea corta de altavoces en la pared del lado de la sala de control
% alt(1,1:8)=1.0182-(0.09:0.18:1.35)*cos(45*pi/180);
% alt(2,1:8)=(0.09:0.18:1.35)*cos(45*pi/180);
% alt(1,9:32)=0;
% alt(2,9:32)=(1.0182)+(0.09:0.18:2*1.44+1.35);
% alt(1,33:40)=(0.09:0.18:1.35)*cos(45*pi/180);
% alt(2,33:40)=(1.0182+3*1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
% alt(1,41:48)=(1.0182)+(0.09:0.18:1.35);
% alt(2,41:48)=(2*1.0182+3*1.44);
% alt(1,49:56)=(1.0182+1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
% alt(2,49:56)=(2*1.0182+3*1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
% alt(1,57:80)=(2*1.0182+1.44);
% alt(2,57:80)=(1.0182)+(2*1.44+1.35:-0.18:0.09);
% alt(1,81:88)=(2*1.0182+1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
% alt(2,81:88)=(1.0182)-(0.09:0.18:1.35)*cos(45*pi/180);
% alt(1,89:96)=(1.0182+1.44)-(0.09:0.18:1.35);
% alt(2,89:96)=0;
nalt=96; % Número de altavoces totales

%FUENTE (la posisicón es arbitraria)
%fte=[3.8;3.2]; %fte 1
%fte=[1.75;7]; %fte 2
%fte=[-2.5;2]; %fte 3
nfte=1;  % Núemero total de fuentes

% Si se quiere representar la distribución de las fuentes y altavoces, hay
% que descomentear la línea siguiente. 
%figure;plot(alt(1,:),alt(2,:),'*'); hold on; plot(fte(1),fte(2),'+'); 


% Datos de los angulos de orientación de cada array lineal (en este caso 6 grupos de 8 altavoces y 2 grupos de 24 altavoces)
% Tomamos como 0 º el que está arriba, sobre el eje X orientado hacia
% valores negativos de X. Lo escribimos en grados.

%tecta=[ones(1,8)*135,ones(1,24)*90,ones(1,8)*45,ones(1,8)*0,ones(1,8)*-45,ones(1,24)*-90,ones(1,8)*-135,ones(1,8)*180];
tecta=[ones(1,8)*45,ones(1,24)*0,-ones(1,8)*45,-ones(1,8)*90,ones(1,8)*-135,ones(1,24)*-180,-ones(1,8)*225,ones(1,8)*90];

fuente=fte;
x = fuente(1);
y = fuente(2);

%Cálculo del angulo relativo fte-altavoz (para determinar si un determinado altavoz debería activarse o no dependiendo de donde esté la fuente que quiere sintenizar)
difX = alt(1,:)-x; % Distancia en la coordenada x
difY = alt(2,:)-y; % % Distancia en la coordenada y
alfa = atan2(difY,difX);      % Mariz con los ángulos de relativos en radianes de cada altavoz con la fuente
alfa=(alfa.*180/pi)+90-tecta; % Matriz con los ángulos anteriores en grados ocn la corrección de la orientación de cada altavoz
[parray,pos]=find(((alfa<90)&(alfa>-90))|((alfa<450)&(alfa>270))); % Buscamos los altavoces que debería estar activos
parray_act=pos;   % esta variable contiene el identificador de cada altavoz de los que estarán activos y de los que hay que calcular las señales que lo alimentan 
%plot(alt(1,pos),alt(2,pos),'o');
%Cálculo de la distancia fte-altavoz
r = sqrt(difX.^2 + difY.^2);                     % Distancia entre la fuente y cada altavoz
%alfa=difX./cos(tecta-pi/2)./r;
%Amplitud
% Para la distancia entre el altavoz y la línea imaginaria interior. Vamos
% a suponer una línea interior con la misma forma que el array. Y que está
% situada a 2/3 del array, contando desde el array hasta el punto central.
% Se tiene que calcular para todo, pero de momento, para una fuente situada
% a la izquierda del array, para los altavoces que están activos (2,3 y 4),
% calculamos la distancia para el centro de la sala y sería: 
% 1.44/2+1.44*cos(45*pi/180)
%r0=ones(1,nalt) *2.88 * (1/2);
%r0=ones(1,nalt) * (1.44/2+1.44*cos(45*pi/180));
r0=ones(1,nalt) *1.28*2;
s0=abs(r);%;.*sin(alfa);
A=sqrt(r0./(r0+s0));
%A=sqrt(1./(s0));
an=A.*cos(alfa*(pi/180))./(sqrt(r)); % En an, se almacena la amplitud con la que habría que corregir la señal de la fuente virtual para emitirla con cada altavoz

%% Retardo
tn=-land.*(fs*(r/c)); % en ts se almacena el retardo que habría que aplicar a la señal de la fuente virtual al emitirla por cada altavoz
t0=round(abs(min(tn)));
tmax=round(abs(max(tn(parray_act))));





%% Filtro salida (en caso de no querer ir retardando y ponderando la señal para cada altavoz, se pueden diseñar filtros que realizan ese efecto).
%% En este caso son filtros FIR de 512 elementos que incorporan el efecto del filtro H_1 calculado anteriormente. 
sn=1;
paux_sign=conv(sn,H_1); %length=length(sn)+length(h)-1
delta=zeros(512-ord,1);
maxDist = max(r);
order = 2^nextpow2(ceil(maxDist/c*fs));
pout_sign_aux=zeros(order,nalt);
for i=parray_act
    paux_out=paux_sign*an(i);
    delta=zeros(512-ord,1);
    delta(round(tn(i)))=1; 
    %pout_sign_aux(:,i)=conv(paux_out,delta);
    pout_sign_aux(round(tn(i)),i)=an(i);
end;

end