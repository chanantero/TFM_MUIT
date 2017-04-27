% Test 1 (25/4/2017). Un solo writer. Comprobamos que todos los par�metros
% de tiempo tienen un comportamiento de acuerdo a la teor�a.

% Crear se�al (Crear el audio file reader, configurarlo y extraer se�al)
readerObj = dsp.AudioFileReader();
readerObj.Filename = 'C:\Users\Rub�n\Music\La Oreja De Van Gogh - El Planeta Imaginario (2016)\La Oreja De Van Gogh - Diciembre.mp3';
Tr_desired = 1.5; % Tiempo de reproducci�n en segundos
readerObj.SamplesPerFrame = ceil(readerObj.SampleRate*Tr_desired);
signal = step(readerObj);
release(readerObj)

% Crear el audio device writer y configurarlo
writerObj = audioDeviceWriter;
writerObj.SampleRate = readerObj.SampleRate;
writerObj.SupportVariableSizeInput = false;
if writerObj.SupportVariableSizeInput
    writerObj.BufferSize = size(signal, 1);
    bufferSize = writerObj.BufferSize;
else
    bufferSize = size(signal, 1);
end
TloadingFrame = size(signal, 1)/writerObj.SampleRate;
TbufferFrame = bufferSize/writerObj.SampleRate;

% Computar
N = 20; % N�mero de iteraciones
t_pausa = 1.5;
t_bp = zeros(N, 1); % Tiempos de incio de procesado
t_ep = zeros(N, 1); % Tiempos de finalizaci�n de procesado
numUnderrun = zeros(N,1); % N�mero de muestras reproducidas en silencio

tref = tic;
for k = 1:N
    t_bp(k) = toc(tref);
    signal = step(readerObj);
    numUnderrun(k) = play(writerObj, signal);
    t_ep(k) = toc(tref);
    pause(t_pausa);
end
t_underrun = numUnderrun/writerObj.SampleRate;

release(writerObj)

% Representar
if writerObj.SupportVariableSizeInput
    numUnderrunBuffers = numUnderrun/writerObj.BufferSize;
else
    numUnderrunBuffers = numUnderrun/size(signal, 1);
end
f = figure;
ax = subplot(3, 2, 1);
plot(t_bp, 'Marker', '.', 'LineWidth', 0.01)
ax.Title.String = 'Time start processing';
ax = subplot(3, 2, 3);
plot(t_ep, 'Marker', '.', 'LineStyle', 'none')
ax.Title.String = 'Time end processing';
ax = subplot(3, 2, 5);
plot(t_ep - t_bp, 'Marker', '.', 'LineStyle', 'none')
ax.Title.String = 'Time processing';
ax = subplot(3, 2, 2);
plot(t_bp(2:end) - t_ep(1:end-1), 'Marker', '.', 'LineStyle', 'none')
ax.NextPlot = 'Add';
plot(diff(t_bp), 'Marker', '.', 'LineStyle', 'none')
ax = subplot(3, 2, 4);
plot(numUnderrunBuffers, 'Marker', '.', 'LineStyle', 'none')
ax.Title.String = 'Number of underrun buffers';

%
[ fig, delayLimits ] = processBufferResults( t_ep, t_underrun, TbufferFrame, TloadingFrame);

% Buffer theory
t_bufferQueueLoad = t_ep; % Tiempos de carga de frame en la cola de buffer
t_lq = TloadingFrame*(1:N)'; % Tiempo de se�al acumulado en la cola de buffer en los tiempos de carga en la cola de buffer (despu�s de la carga)

t_delay = 0;
t_br = t_delay + cumsum([t_underrun(1); t_underrun(2:end) + TloadingFrame]); % Tiempos de inicio de reproducci�n
t_er = t_br + TloadingFrame; % Tiempos de fin de reproducci�n

t_frameChange = (t_br(1):TloadingFrame:t_er(end))'; % Tiempos de cambio de frame (carga y lectura de buffer)
t_lr = interp1(t_er, (1:N)*TbufferFrame, t_frameChange, 'previous', 0); % Tiempo de se�al reproducido hasta el momento en cada cambio de buffer

% Tiempo que pasa entre un cambio de de frame y el momento en que ya hay
% suficiente se�al en la cola de buffer como para reproducir otro
t_enough = interp1([0; t_lq], [0; t_bufferQueueLoad], t_lr + TbufferFrame, 'next', Inf); % Tiempos en los que ya hay suficiente se�al en la cola de buffer como para reproducir otro frame
t_waiting = t_enough - t_frameChange;

% Condici�n 1: La reproducci�n de un frame no vac�o siempre ha de producirse cuando hay
% suficiente se�al en la cola de buffer como para reproducirse. Es decir,
% t_waiting ha de ser negativo siempre que una reproducci�n se inicie
% 
% Condici�n 2: La reproducci�n de un frame vac�o siempre ha de producirse cuando todav�a
% no hay suficiente se�al en la cola de buffer como para reproducirse. Es
% decir, t_waiting ha de ser positivo siempre que una reproducci�n de frame
% vac�o (underrun) ocurra
rbFlag = ismember(t_frameChange, t_br); % �ndices de cambios de frame en los que se inicia una reproducci�n
underrunFlag = ~rbFlag; % �ndices de cambios de frame en los que se inicia un underrun (frame vac�o)

lowLimitDelay = max(t_waiting(rbFlag)); % L�mite inferior de retardo que puede a�adirse para que se cumpla la condici�n de causalidad (condici�n 1)
upLimitDelay = min(t_waiting(underrunFlag)); % L�mite superior de retardo que puede a�adirse para que se cumpla la condici�n 2
assert(upLimitDelay >= lowLimitDelay, 'Error!!!')
differ = upLimitDelay - lowLimitDelay; % Margen de indeterminaci�n

taddDelay = max(0, lowLimitDelay); % Retardo que hay que a�adirse al tiempo de inicio de reproducci�n actual
tSubstractDelay = max(0, -upLimitDelay); % Retardo que hay que sustraer al tiempo de inicio de reproducci�n actual

% Representar
% La se�al acumulada en la cola de buffer evoluciona discontinuamente. Hay
% que representarlo a escalones.
x = kron(t_bufferQueueLoad, [1; 1]);
base = 0; % Cantidad de se�al que hab�a antes de la primera carga de cola
y = reshape([[base; t_lq(1:end-1)], t_lq]', 2*N, 1);

plot(t_frameChange, t_lr, x, y)


%%
% Test. �C�mo se comporta la reproducci�n en un deviceWriter cuando
% introduzco pause despu�s de volcar la se�al al buffer de escritura? �Se
% espera el tiempo que indico en pause para reproducir? �Y si no se espera,
% qu� ocurre con el numUnderrun devuelto?
% Cosas descubiertas:
% - El pause funciona tal y como se supone, no hay problema con eso.
% - Cuando SupportVariableSizeInput es false, numUnderrun es siempre un
% m�ltiplo del n�mero de muestras de entrada al writer object
% - Cuando SupportVariableSizeInput es true, el numUnderrun se reinicia a 0
% cuando por primera vez cambio el tama�o del input
% - Cuando SupportVariableSizeInput es true, numUnderrun es siempre un
% m�ltiplo de BufferSize

% Crear el audio device writer y configurarlo
writerObj = audioDeviceWriter;
writerObj.SupportVariableSizeInput = true;
writerObj.BufferSize = 36000;

% Crear el audio device reader y configurarlo
readerObj = dsp.AudioFileReader();
readerObj2 = dsp.AudioFileReader();
readerObj.Filename = 'C:\Users\Rub�n\Music\La Oreja De Van Gogh - El Planeta Imaginario (2016)\La Oreja De Van Gogh - Diciembre.mp3';
readerObj2.Filename = 'C:\Users\Rub�n\Music\Varias\Glee Defying gravity season 5 lyrics (128  kbps).mp3';
readerObj.SamplesPerFrame = readerObj.SampleRate;
readerObj2.SamplesPerFrame = round(readerObj2.SampleRate/4);

signal1 = step(readerObj);
signal2 = step(readerObj2);
% signal = zeros(readerObj.SampleRate, 2);
numUnderrun = zeros(3,1);
buffer = 0;
play(writerObj, zeros(2,2));
play(writerObj, zeros(1, 2));

numUnderrun(1) = play(writerObj, signal1);
buffer = rem(buffer + readerObj.SampleRate, writerObj.BufferSize);

numUnderrun(2) = play(writerObj, signal1);
buffer = rem(buffer + readerObj.SampleRate, writerObj.BufferSize);

pause(4)
numUnderrun(3) = play(writerObj, signal2);
buffer = rem(buffer + readerObj2.SampleRate, writerObj.BufferSize);



double(numUnderrun)/writerObj.SampleRate


release(writerObj)
release(readerObj)