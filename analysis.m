%% An�lisis de se�ales de micr�fonos

s_micro; % Se�ales de micr�fonos. Cada columna es un canal

% En teor�a, las se�ales son sinusoidales, por lo que, sabiendo la
% frecuencia, solo nos queda una inc�gnita que encontrar: la envolvente
% Tambi�n conocemos la envolvente de las fuentes
iq = real2IQ(s_micro, SampleRate, freq);

x = cos(