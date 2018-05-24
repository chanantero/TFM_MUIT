%% Basic simulation
% Reduced version of the simulations

%% System parameters

% Constants
c = 340; % Sound velocity (m/s)
fs = 20000; % Sample frequency (samples/s)

% Noise source coefficient
amplitude = 1;
phase = 0;

% % Creation of frequency filter
% magnFiltOrder = 2^(12);
% hilbertFiltOrder = 2^(12);
% [freqFilter, freqFiltDelay] = getFrequencyFilter( magnFiltOrder, hilbertFiltOrder, fs );
% freqFilterLength = length(freqFilter);

% % Position of loudspeakers of the WFS array
% d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
% nb = 8; % Bottom and upper sides of the octogon (2 sides)
% nd = 8; % Diagonal sides of the octogon (4 sides)
% nl = 24; % Lateral side of the octogon (2 sides)
% betabd = 45; % Deviation angle between bottom/upper and diagonal sides
% [ x, y, alfa ] = octogon(d, nb, nd, nl, betabd);
% z = zeros(numel(x), 1);
% WFSpositions = [x, y, z];
% loudspeakersOrientation = [cosd(alfa), sind(alfa), zeros(numel(alfa), 1)];
% numWFS = size(WFSpositions, 1); % 96  
% centerX = (max(x) + min(x))/2;
% centerY = (max(y) + min(y))/2;

% Microphone position
% recPosition = [centerX, centerY, 0];

% Positions of the noise source
% NSposition = [centerX + 5, centerY, 0];    
NSposition = [0.1 2.3 1.65];

durSign = 1; % Duration of signal
t = (0:ceil(durSign*fs)-1)/fs;
NSsignal = chirp(t, 20, durSign, 940);

%% Room characteristics and impulse response of chamber
[h_array,alt,tecta,L,mallado_x,mallado_y,h_ad_sources] = SalaGtac(1.65,0.4,0.4,0,1000,fs,c,fte);
WFSpositions = [alt' 1.65*ones(96, 1)];
roomDim = L;
[X, Y] = ndgrid(mallado_x, mallado_y);
recPosition = [X(:), Y(:), 1.65*ones(numel(X), 1)];
loudspeakersOrientation = [cosd(tecta'-90), sind(tecta' - 90), zeros(96, 1)];
numRec = size(recPosition, 1);

numReverbTime = 1;
beta = 0; % Average reflection coefficient of the walls of the chamber
Beta = beta(:) * [1 1 1 1 1 1];
r = recPosition; % Receiver position [x y z] (m)
wfsPos = WFSpositions;
nsPos = NSposition;
maxX = max([nsPos(:, 1); wfsPos(:, 1)]);
maxY = max([nsPos(:, 2); wfsPos(:, 2)]);

% Adjust the number of samples of impulse responses
numSampIR = 1000; % Number of samples

WFSfilterLength = numSampIR*2;

WFS_IR = zeros(numSampIR, numWFS, numRec);
for k = 1:numWFS
    WFS_IR(:, k, :) = permute(rir_generator(c, fs, r, wfsPos(k, :), roomDim, Beta(rt, :), numSampIR), [2 3 1]);
end

NS_IR = rir_generator(c, fs, r, nsPos, roomDim, Beta, numSampIR)';

ind = 80;
[indX, indY] = ind2sub([length(mallado_x), length(mallado_y)], ind);
a = h_array(:,:,indX,indY);
b = WFS_IR(:, :, ind);
isequal(a,b)

isequal(NS_IR(:, ind), h_ad_sources(:, 1, indX, indY))

%% Simulation
% preN = freqFiltDelay;
% postN = numSampIR - 1 + freqFiltDelay - 1;
% x = [zeros(1, preN), NSsignal, zeros(1, postN)];
% numSamp = length(x);

% Generate signals of WFS array loudspeakers

% Calculate delay and attenuation
relPos = WFSpositions - repmat(NSposition, numWFS, 1);
distances = sqrt(sum(relPos.^2, 2));
delays = distances/c;
r0 = 1.44*(0.5+cosd(45));
A = sqrt(r0./(r0 + distances));
r0_Miguel = 1.28*2;
A_Miguel=sqrt(r0_Miguel./(r0_Miguel+distances));
cosAlfa = dot(relPos, loudspeakersOrientation, 2)./distances;
attenuations = -A_Miguel.*cosAlfa./sqrt(distances)*d;
attenuations(cosAlfa < 0) = 0;

[filtros_array,an,tn,activo_array]=WFS_DrivingSignals(alt,fte(fuente_ruido,:),c,fs);

% tecta=[ones(1,8)*45,ones(1,24)*0,-ones(1,8)*45,-ones(1,8)*90,ones(1,8)*-135,ones(1,24)*-180,-ones(1,8)*225,ones(1,8)*90];
% x = NSposition(1);
% y = NSposition(2);
% %Cálculo del angulo relativo fte-altavoz (para determinar si un determinado altavoz debería activarse o no dependiendo de donde esté la fuente que quiere sintenizar)
% difX = alt(1,:)-x; % Distancia en la coordenada x
% difY = alt(2,:)-y; % % Distancia en la coordenada y
% alfa = atan2(difY,difX);      % Mariz con los ángulos de relativos en radianes de cada altavoz con la fuente
% alfa=(alfa.*180/pi)+90-tecta; % Matriz con los ángulos anteriores en grados ocn la corrección de la orientación de cada altavoz
% [parray,pos]=find(((alfa<90)&(alfa>-90))|((alfa<450)&(alfa>270))); % Buscamos los altavoces que debería estar activos
% parray_act=pos;   % esta variable contiene el identificador de cada altavoz de los que estarán activos y de los que hay que calcular las señales que lo alimentan 
% %plot(alt(1,pos),alt(2,pos),'o');
% %Cálculo de la distancia fte-altavoz
% r = sqrt(difX.^2 + difY.^2);
% nalt=96; % Número de altavoces totales
% r0_Miguel =ones(1,nalt) *1.28*2;
% s0=abs(r);%;.*sin(alfa);
% A_Miguel=sqrt(r0_Miguel./(r0_Miguel+s0));
% an_Miguel =A_Miguel.*cos(alfa*(pi/180))./(sqrt(r)); % En an, se almacena la amplitud con la que habría que corregir la señal de la fuente virtual para emitirla con cada altavoz

% an_Ruben = A'.*cos(alfa*(pi/180))./(sqrt(r));
% attenuations = -A.*cosAlfa./sqrt(distances);
% [an_Ruben', attenuations]

% Principales diferencias entre el cálculo de unas atenuaciones y otras.
% - Miguel no usa el escalado por la separación entre altavoces.
% - Tampoco usa el mismo r0
% - No multiplica por -1

% Los retardos los devuelve en muestras en vez de en segundos
% [tn'/fs, delays]

% Definimoa la señal de ruido
durSign = 1; % Duration of signal
t = (0:ceil(durSign*fs)-1)/fs;
f = 200;
NSsignal = sin(2*pi*f*t)';
x = NSsignal;
% x=sin(2*pi*0.06*[0:1999])';
numSamp = length(x);
% Definimos la señal de la fuente virtual;
xv=-0.0957*x;
%xv=-x;
xv(1:999)=0; % empieza a actuar en la muestra 1000;

[POT_ad,POT_ar]=generamapa(xv,h_array,filtros_array,activo_array,x,h_ad_sources(:,fuente_ruido,:,:));

% Representamos el mapa de presiones
 dv=50; % Duración de la ventana para promediar al potencia
 despv=25;  % Desplazamiento del enventanado para realizar el calculo de la
% potencia

% dibujapot(POT_ad+POT_ar ,L, alt, fte(fuente_ruido,:), mallado_x, mallado_y,dv,despv);


% Create basic filter impulse response (delta)
indDelta = floor(delays*fs) + 1;
NsFilt = max(floor(delays(:)*fs)) + 1 + freqFilterLength;
    
filtersWFS = zeros(numWFS, NsFilt);
for ss = 1:numWFS
    if attenuations(ss) ~= 0
        filtersWFS(ss, indDelta(ss)) = attenuations(ss);
    end
end

% % Convolute it with the frequency filter
% filtersWFS = fftfilt(freqFilter, filtersWFS')';

% Apply it to noise source signal
wfsSignals = fftfilt(filtersWFS', x')';
% % Compensate delay
% wfsSignals = [wfsSignals(:, freqFiltDelay + 1:end), zeros(numWFS, freqFiltDelay)];

% Simulate
acousticPaths = [WFS_IR, permute(NS_IR, [1 3 2])];
sourceSignals = [wfsSignals', x];

recSignals = zeros(numSamp, numRec);
for rec = 1:numRec
    recSignals(:, rec) = sum(fftfilt(acousticPaths(:, :, rec), sourceSignals), 2);
end

% Simulate with only the noise source
recSignalsNS = zeros(numSamp, numRec);
for rec = 1:numRec
    recSignalsNS(:, rec) = sum(fftfilt(acousticPaths(:, end, rec), sourceSignals(:,end)), 2);
end

%% Visualization
t = (0:numSamp-1)/fs;

recSignalsWFS = recSignals - recSignalsNS;
recSignals2 = recSignalsNS + recSignalsWFS*sqrt(f/c);
plot(t, recSignalsNS(:, 1), t, recSignals2(:,1), t, recSignalsWFS(:,1)*sqrt(f/c))

plot(t, recSignalsNS(:, 1), t, recSignals(:,1), t, - recSignalsNS(:, 1) + recSignals(:, 1))
plot(t, recSignalsNS(:, 1), t, POT_ad(:, 1))
plot(t, recSignals(:, 1), t, POT_ar(:,1) + POT_ad(:,1))
plot(t, recSignalsWFS(:, 1), t, POT_ar(:, 1))
plot(t, POT_ad(:, 1), t, POT_ar(:,1), t, 0.98*POT_ar(:,1) + POT_ad(:,1))



