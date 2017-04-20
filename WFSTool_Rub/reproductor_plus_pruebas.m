player = reproductor_plus();

comMat = [1 1];
player.setProps('comMatrix', comMat);
player.setProps('audioFileName', 'C:\Users\Rubén\Music\Ariana Grande - Dangerous Woman\01 - Moonlight.mp3', 1);
player.setProps('getDelayFun', @() [0; 0], [1 1]);
player.setProps('getAttenFun', @() [1; 1], [1 1]);
player.setProps('getDelayFun', @() [0; 0], [1 2]);
player.setProps('getAttenFun', @() [1; 1], [1 2]);
player.setProps('frameDuration', 2);

player.setProps('device', 'Altavoces (Dispositivo de High Definition Audio)', 1);
player.setProps('device', 'Controlador primario de sonido', 2);

setup(player, [], []);
delays = player.testDelayBetweenDevices();

order.action = 'stop';
player.executeOrder(order);

order.action = 'play';
player.executeOrder(order);

% Test. ¿Cómo se comporta la reproducción en un deviceWriter cuando
% introduzco pause después de volcar la señal al buffer de escritura? ¿Se
% espera el tiempo que indico en pause para reproducir? ¿Y si no se espera,
% qué ocurre con el numUnderrun devuelto?
% Cosas descubiertas:
% - El pause funciona tal y como se supone, no hay problema con eso.
% - Cuando SupportVariableSizeInput es false, numUnderrun es siempre un
% múltiplo del número de muestras de entrada al writer object
% - Cuando SupportVariableSizeInput es true, el numUnderrun se reinicia a 0
% cuando por primera vez cambio el tamaño del inputç
% - Cuando SupportVariableSizeInput es true, numUnderrun es siempre un
% múltiplo de BufferSize
writerObj = audioDeviceWriter;
writerObj.SupportVariableSizeInput = true;
readerObj = dsp.AudioFileReader();
readerObj.Filename = 'C:\Users\Rubén\Music\Ariana Grande - Dangerous Woman\01 - Moonlight.mp3';
readerObj.SamplesPerFrame = readerObj.SampleRate;

% signal = step(readerObj);
% signal = zeros(readerObj.SampleRate, 2);
numUnderrun = zeros(3,1);
play(writerObj, zeros(2,2));
play(writerObj, zeros(1, 2));
numUnderrun(1) = play(writerObj, signal);
numUnderrun(2) = play(writerObj, signal);
pause(4)
numUnderrun(3) = play(writerObj, zeros(1, 2));

double(numUnderrun)/writerObj.SampleRate


release(writerObj)
release(readerObj)