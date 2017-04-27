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
t_pausa = 2;
t_bp = zeros(N, 1); % Tiempos de incio de procesado
t_ep = zeros(N, 1); % Tiempos de finalizaci�n de procesado
numUnderrun = zeros(N,1); % N�mero de muestras reproducidas en silencio

tref = tic;
for k = 1:N
    t_bp(k) = toc(tref);
    numUnderrun(k) = play(writerObj, signal);
    t_ep(k) = toc(tref);
    pause(t_pausa);
end

release(writerObj)

% Procesar resultados
t_underrun = numUnderrun/writerObj.SampleRate;
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


% Buffer theory
t_bufferQueueLoad = t_ep; % Tiempos de carga de frame en la cola de buffer
t_lq = TloadingFrame*(1:N)'; % Tiempo de se�al acumulado en la cola de buffer en los tiempos de carga en la cola de buffer

t_delay = 0;
t_br = t_delay + cumsum([t_underrun(1); t_underrun(2:end) + TloadingFrame]); % Tiempos de inicio de reproducci�n
t_er = t_br + TloadingFrame; % Tiempos de fin de reproducci�n

t_frameChange = (t_br(1):TloadingFrame:t_er(end))'; % Tiempos de cambio de buffer (carga y lectura de buffer)
t_lr = zeros(size(t_frameChange)); % Tiempo de se�al reproducido hasta el momento en cada cambio de buffer
for k = 1:numel(t_frameChange)
    ind = find(t_frameChange(k) >= t_er, 1, 'last');
    if isempty(ind)
        ind = 0;
    end
    t_lr(k) = ind*TloadingFrame;
end


% Buffer queue
tl_q_x = kron(t_bufferQueueLoad, [1; 1]);
tl_q_y = reshape([[0; t_lq(1:end-1)], t_lq]', 2*N, 1); % cumsum(reshape([zeros(1, N); Tr*ones(1, N)], 2*N, 1));

t_r_x = t_frameChange;
t_r_y = t_lr;

plot(t_r_x, t_r_y, tl_q_x, tl_q_y)

% Condiciones del tiempo de retardo entre la finalizaci�n del procesado y
% el inicio de la reproducci�n

% Tiempo de inicio de reproducci�n es siempre igual o posterior al tiempo
% de carga en la cola de buffer
delayRep = t_br - t_ep; % Tiempo que pasa entre la carga en la cola de buffer y la reproducci�n de ese frame
taddDelay = max(0, -min(delayRep)); % Retardo que hay que a�adirse al tiempo de inicio de reproducci�n actual

% No deber� haber underrun si hay muestras disponibles en la cola de buffer
% Tiempos de carga de buffer vac�o (underrun)
iniRep_log = ismember(t_frameChange, t_br); % �ndices de tiempos de inicio de reproducci�n
iniEmptyBuffer_log = true(numel(t_frameChange), 1); % �ndices de tiempos de inicio de underrun
iniEmptyBuffer_log(iniRep_log) = false;
t_b_EmptyBuffer = t_frameChange(iniEmptyBuffer_log);
acumTime_reprod = t_lr(iniEmptyBuffer_log);

% Cu�nta se�al acumulada hay en la cola en esos tiempos?
acumTime_queue = zeros(numel(t_b_EmptyBuffer), 1);
for k = 1:numel(t_b_EmptyBuffer)
    ind = find(t_b_EmptyBuffer(k) >= t_bufferQueueLoad, 1, 'last');
    if isempty(ind)
        acumTime_queue(k) = 0;
    else
        acumTime_queue(k) = t_lq(ind);
    end
end

acumTime_queue - acumTime_reprod < bufferSize/writerObj.SampleRate; % Condici�n que ha de cumplirse

% Cu�l es la reducci�n m�nima de retardo necesaria para cumplir con la
% condici�n?
% Tiempo que pasa entre el cambio de frame de reproducci�n y el tiempo en
% que la cola de buffer tiene suficiente se�al como para seguir
% reproduciendo
for k = 1:numel(t_frameChange)
    ind = find(t_lq >= t_lr(k) + TbufferFrame, 1, 'first');
    a(k) = t_bufferQueueLoad(ind) - t_frameChange(k);
end


t_enough = zeros(numel(acumTime_reprod), 1); % Tiempo que hace que ya hay suficiente se�al en la cola de buffer para reproducir
for k = 1:numel(acumTime_reprod)
    ind = find(t_lq <= acumTime_reprod(k), 1, 'last');
    t_enough(k) = t_lq(ind) - acumTime_reprod(k);
end
tsubstractDelay = max(t_enough);



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