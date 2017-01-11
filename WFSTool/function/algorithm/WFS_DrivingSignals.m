function [parray_act,an,tn]=WFS_DrivingSignals(vsx,vsy)
%%%% The aim of this code is to compute amplitudes and delays that correspond
%%%% to every loduspeaker. The loudspeakers belong to a 96 array with
%%%% octogonal shape.
%%%% The loduspeaker positions are defined in 2x96 matrix "alt". One row
%%%% for each cartesian coordinate and one column for each loudspeaker. 
%%%% The variable fte gives the virtual position (it contains x and y source coordinates). 
%%%% Morover, code computes by using trigonometry the active loduspeakers. 
%%%% The selection of those active loduspeakers depend on the angle that the 
%%%% sound coming from the virtual source had upon entry. 
%%%% Output variables:
%%%% parray_act -> identifies the active loduspeakers. Example:
%%%% if all of them are active parra_act=[1:96];
%%%% an-> vector with the gain of the virtual source signal at every loudspeaker.
%%%% tn-> vector with the delays of the the virtual source signal at every loudspeaker.
%%%% Finally, pout_sign_aux is a matrix with 512 rows and so many columns
%%%% as loudspeakers has the array. It stores the coefficients of a 512-taps FIR filter
%%%% for each loudspeaker using the gain and delays previously computed. 
%%%% This way we can use those filters to filter the virtual source signal 
%%%% before fed the loudspeakers
%%%% The code by itself does not work with signals. That means it does not synthesize
%%%% the sound field. It only computes the correction factors (gains and delays) 
%%%% required to generate the signal at every loudspeaker. Given a signal, 
%%%% those factors allow to obtain the signals that fed every loudpseaker. 
%%%% If the designed filters are used, we only had to filter the signal for every filter before
%%%% to fed the corresponding loudspeaker.  

c=343; % Sound speed
comp=sqrt(-1); % j
ang_max=90*(pi/180); % Angulo máximo de incidencia del sonido (90? en radianes
ord=60;
Ax=0.18; % Separación entre altavoces
FS=44100; % Frecuencia de trabajo ²ÉÑùÂÊ
land=-1; % porque est?fuera del array de altavoces
ali_frec=c/(2*max(Ax)*sin(ang_max)); % frecuencia de alias
ali_frec=ceil(ali_frec); % frecuencia de alias redondeada al alza

% Diseño de filtro para la reconstrucción de la señal en campo cercano
% WFS podría funcionar sin este filtro
frec_filtro=linspace(0,ceil(FS/2)+1,(ceil(FS/2)/40)+1);
longit=length(frec_filtro)-mod(length(frec_filtro),2);
amp_filtro=sqrt(-(frec_filtro(1:ceil((ali_frec+1000)/40)))/(c*comp));
vector(1:ceil((ali_frec+1000)/40))=amp_filtro;
vector(ceil((ali_frec+1000)/40)+1:longit)=amp_filtro(end);
H_1=firls(ord,frec_filtro(1:longit)/frec_filtro(longit),abs(vector)); % H_1 contiene la respuesta al impulso del filtro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Definición de la geometría del array (que se corresponde con el de la sala GTAC)
%% Coordenadas x,y de las posiciones de los altavoces (el eje x sería paralelo a la pared larga) 
%% x=0 se corresponde con la línea larga de los altavoces en la pared de enfrente de la puerta
%% y=0 se corresponde con la línea corta de altavoces en la pared del lado de la sala de control
[alt, tecta] = generate_array_PosOrient();
nalt = size(alt, 2); %96; % number of loudspeakers in the array

%FUENTE (la posisicón es arbitraria)
%fte=[3.8;3.2]; %fte 1
%fte=[1.75;7]; %fte 2
fte=[vsx;vsy]; %fte 3
nfte=1;  % number of fte

% Si se quiere representar la distribución de las fuentes y altavoces, hay
% que descomentear la línea siguiente. 
% figure;plot(alt(1,:),alt(2,:),'*'); hold on; plot(fte(1),fte(2),'+'); 

% Datos de los angulos de orientación de cada array lineal (en este caso 6 grupos de 8 altavoces y 2 grupos de 24 altavoces)
% Tomamos como 0 ?el que est?arriba, sobre el eje X orientado hacia
% valores negativos de X. Lo escribimos en grados.

fuente=fte;
x = fuente(1);
y = fuente(2);

%Cálculo del angulo relativo fte-altavoz (para determinar si un determinado altavoz debería activarse o no dependiendo de donde est?la fuente que quiere sintenizar)
% Un altavoz se activa si la línea que une la fuente virtual con el altavoz y
% la dirección principal en la que el altavoz emite ("broadside") forman un
% ángulo menor a 90º
% Particular para el caso 2D
difX = alt(1,:)-x; % Distancia en la coordenada x
difY = alt(2,:)-y; % % Distancia en la coordenada y
alfa = atan2d(difY,difX) - tecta; % Matriz con los ángulos anteriores en grados con la corrección de la orientación de cada altavoz 
parray_act = find(((alfa<90)&(alfa>-90))|((alfa<450)&(alfa>270))); % Buscamos los altavoces que debería estar activos. % esta variable contiene el identificador de cada altavoz de los que estarán activos y de los que hay que calcular las señales que lo alimentan

% Caso 3D
[coord3D, orient] = generate_array_3D(); % Posición y orientación de altavoces
fuente3D = [x, y, 0];
relPos = coord3D - repmat(fuente3D, nalt, 1); % Posición relativa respecto a la fuente
r = sqrt(sum(relPos.^2, 2));
normOrient = sqrt(sum(orient.^2, 2)); % Ideally 1
cosAlfa = dot(relPos, orient, 2)./(r.*normOrient); % Coseno del ángulo entre la línea que une la fuente virtual con cada altavoz y
% la dirección principal en la que el altavoz emite ("broadside")
parray_act3D = find(cosAlfa > 0)'; % Buscamos los altavoces que debería estar activos. % esta variable contiene el identificador de cada altavoz de los que estarán activos y de los que hay que calcular las señales que lo alimentan

%Cálculo de la distancia fte-altavoz
r = sqrt(difX.^2 + difY.^2);  % Distancia entre la fuente y cada altavoz
%alfa=difX./cos(tecta-pi/2)./r;
%Amplitud
% Para la distancia entre el altavoz y la línea imaginaria interior. Vamos
% a suponer una línea interior con la misma forma que el array. Y que est?
% situada a 2/3 del array, contando desde el array hasta el punto central.
% Se tiene que calcular para todo, pero de momento, para una fuente situada
% a la izquierda del array, para los altavoces que están activos (2,3 y 4),
% calculamos la distancia para el centro de la sala y sería: 
% 1.44/2+1.44*cos(45*pi/180)
%r0=ones(1,nalt) *2.88 * (1/2);
r0=ones(1,nalt) * (1.44/2+1.44*cos(45*pi/180));
s0=abs(r);%.*sin(alfa);
A=sqrt(r0./(r0+s0));
an=A.*cos(alfa*(pi/180))./(sqrt(r)); % En an, se almacena la amplitud con la que habría que corregir la señal de la fuente virtual para emitirla con cada altavoz

%% Retardo
tn=-land.*(FS*(r/c)); % en ts se almacena el retardo que habría que aplicar a la señal de la fuente virtual al emitirla por cada altavoz
t0=round(abs(min(tn)));
tmax=round(abs(max(tn(parray_act))));

%% Filtro salida (en caso de no querer ir retardando y ponderando la señal para cada altavoz, se pueden diseñar filtros que realizan ese efecto).
%% En este caso son filtros FIR de 512 elementos que incorporan el efecto del filtro H_1 calculado anteriormente. 
% sn=1;
% paux_sign=conv(sn,H_1); %length=length(sn)+length(h)-1
% delta=zeros(512-ord,1);
% pout_sign_aux=zeros(512,nalt);
% for i = parray_act
%     paux_out=paux_sign*an(i);
%     delta=zeros(512-ord,1);
%     delta(round(tn(i)))=1; 
%     pout_sign_aux(:,i)=conv(paux_out,delta);
% end;
