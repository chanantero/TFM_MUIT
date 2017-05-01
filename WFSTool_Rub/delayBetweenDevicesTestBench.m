% Test 1 (25/4/2017). Un solo writer. Comprobamos que todos los parámetros
% de tiempo tienen un comportamiento de acuerdo a la teoría.

% Crear señal (Crear el audio file reader, configurarlo y extraer señal)
readerObj = dsp.AudioFileReader();
readerObj.Filename = 'C:\Users\Rubén\Music\La Oreja De Van Gogh - El Planeta Imaginario (2016)\La Oreja De Van Gogh - Diciembre.mp3';
Tr_desired = 2; % Tiempo de reproducción en segundos
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
N = 40; % Número de iteraciones
t_pausa = 1.9; increase = 0.1;
t_bp = zeros(N, 1); % Tiempos de incio de procesado
t_ep = zeros(N, 1); % Tiempos de finalización de procesado
numUnderrun = zeros(N,1); % Número de muestras reproducidas en silencio

tref = tic;
for k = 1:N
    t_bp(k) = toc(tref);
    signal = step(readerObj);
    numUnderrun(k) = play(writerObj, signal);
    if numUnderrun(k) == 0
        t_pausa = t_pausa + increase;
    elseif numUnderrun(k) >= 2*bufferSize
        t_pausa = t_pausa - increase;
    end
    t_ep(k) = toc(tref);
    pause(t_pausa);
end
t_underrun = numUnderrun/writerObj.SampleRate;
numUnderrunBuffers = numUnderrun/bufferSize;

release(writerObj)

% Representar
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

[ fig, delayLimits ] = processBufferResults( t_ep, numUnderrunBuffers, TloadingFrame, TbufferFrame);



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