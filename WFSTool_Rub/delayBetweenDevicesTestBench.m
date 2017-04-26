% Test 1 (25/4/2017). Un solo writer. Comprobamos que todos los parámetros
% de tiempo tienen un comportamiento de acuerdo a la teoría.

% Crear señal (Crear el audio file reader, configurarlo y extraer señal)
readerObj = dsp.AudioFileReader();
readerObj.Filename = 'C:\Users\Rubén\Music\La Oreja De Van Gogh - El Planeta Imaginario (2016)\La Oreja De Van Gogh - Diciembre.mp3';
Tr_desired = 1.5; % Tiempo de reproducción en segundos
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
Tr = size(signal, 1)/writerObj.SampleRate;

% Computar
N = 20; % Número de iteraciones
t_pausa = 2;
t_bp = zeros(N, 1); % Tiempos de incio de procesado
t_ep = zeros(N, 1); % Tiempos de finalización de procesado
numUnderrun = zeros(N,1); % Número de muestras reproducidas en silencio

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
t_delay = 0;
t_br = t_delay + cumsum([t_underrun(1); t_underrun(2:end) + Tr]); % Tiempos de inicio de reproducción
t_er = t_br + Tr; % Tiempos de fin de reproducción

t_bufferLoad = t_br(1):Tr:t_br(end); % Tiempos de inicio de carga y lectura de buffer
t_lr = zeros(size(t_bufferLoad)); % Tiempo de señal reproducido hasta el momento en cada inicio de carga y lectura de buffer
for k = 1:numel(t_bufferLoad)
    ind = find(t_bufferLoad >= t_er, 1, 'last');
    if isempty(ind)
        ind = 0;
    end
    t_lr(k) = ind*Tr;
end


% Buffer queue
tl_q_x = kron(t_ep, [1; 1]);
tl_q_y = cumsum(reshape([zeros(1, N); Tr*ones(1, N)], 2*N, 1));

t_r_x = reshape([t_br, t_er]', 2*N, 1);
t_r_y = cumsum(reshape([zeros(1, N); Tr*ones(1, N)], 2*N, 1));

plot(t_r_x, t_r_y, tl_q_x, tl_q_y)

% Condiciones del tiempo de retardo entre la finalización del procesado y
% el inicio de la reproducción

% Tiempo de inicio de reproducción es siempre igual o posterior al tiempo
% de finalización del procesado (t_ep <= tf_r)
taddDelay = max(0, max(t_ep - t_br)); % Retardo que hay que añadirse al tiempo de inicio de reproducción actual

% Si asumimos que el tiempo de procesado es igual al tiempo en el que las
% muestras están disponibles en la cola de buffer, no deberá haber tiempo
% tiempo underrun si hay muestras disponibles
flag_Underrun = t_underrun(1:end) > 0;

bufferLoad = t_br(1):Tr:t_br(end); % Tiempos de carga de buffer
% Tiempos de carga de buffer vacío (underrun)
iniRep_log = ismember(bufferLoad, t_br); % Índices de tiempos de inicio de reproducción
iniEmptyBuffer_log = true(numel(bufferLoad), 1); % Índices de tiempos de inicio de reproducción
iniEmptyBuffer_log(iniRep_log) = false;
tb_EmptyBuffer = bufferLoad(iniEmptyBuffer_log);

% Cuánta señal acumulada hay en esos tiempos?
ind = zeros(numel(tb_EmptyBuffer), 1);
for k = 1:numel(tb_EmptyBuffer)
    ind(k) = find(tb_EmptyBuffer(k) >= t_ep, 1, 'last');
end
acumTime_reprod = (find(flag_Underrun) - 1)*Tr;
acumTime_queue = ind*Tr;
acumTime - acumTime_reprod < bufferSize/writerObj.SampleRate; % Condición que ha de cumplirse
% Cuál es la reducción mínima de retardo necesaria para cumplir con la
% condición?

%%
% Test. ¿Cómo se comporta la reproducción en un deviceWriter cuando
% introduzco pause después de volcar la señal al buffer de escritura? ¿Se
% espera el tiempo que indico en pause para reproducir? ¿Y si no se espera,
% qué ocurre con el numUnderrun devuelto?
% Cosas descubiertas:
% - El pause funciona tal y como se supone, no hay problema con eso.
% - Cuando SupportVariableSizeInput es false, numUnderrun es siempre un
% múltiplo del número de muestras de entrada al writer object
% - Cuando SupportVariableSizeInput es true, el numUnderrun se reinicia a 0
% cuando por primera vez cambio el tamaño del input
% - Cuando SupportVariableSizeInput es true, numUnderrun es siempre un
% múltiplo de BufferSize

% Crear el audio device writer y configurarlo
writerObj = audioDeviceWriter;
writerObj.SupportVariableSizeInput = true;
writerObj.BufferSize = 36000;

% Crear el audio device reader y configurarlo
readerObj = dsp.AudioFileReader();
readerObj2 = dsp.AudioFileReader();
readerObj.Filename = 'C:\Users\Rubén\Music\La Oreja De Van Gogh - El Planeta Imaginario (2016)\La Oreja De Van Gogh - Diciembre.mp3';
readerObj2.Filename = 'C:\Users\Rubén\Music\Varias\Glee Defying gravity season 5 lyrics (128  kbps).mp3';
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