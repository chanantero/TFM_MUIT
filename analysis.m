%% Análisis de señales de micrófonos

s_micro; % Señales de micrófonos. Cada columna es un canal

% En teoría, las señales son sinusoidales, por lo que, sabiendo la
% frecuencia, solo nos queda una incógnita que encontrar: la envolvente
% También conocemos la envolvente de las fuentes
iq = real2IQ(s_micro, SampleRate, freq);

x = cos(