function [ fig, delayLimits ] = processBufferResultsInvariableInput( t_bufferQueueLoad, t_underrun, T_frame)
% Invariable input size

N = numel(t_bufferQueueLoad); % Número de frames cargados. numel(t_bufferQueueLoad) == numel(t_underrun)

% Buffer theory
% t_bufferQueueLoad: Tiempos de carga de frame en la cola de buffer
t_lq = T_frame*(1:N)'; % Tiempo de señal acumulado en la cola de buffer en los tiempos de carga en la cola de buffer (después de la carga)

t_delay = 0;
t_br = t_delay + cumsum([t_underrun(1); t_underrun(2:end) + T_frame]); % Tiempos de inicio de reproducción
t_er = t_br + T_frame; % Tiempos de fin de reproducción

t_frameChange = (t_br(1):T_frame:t_er(end))'; % Tiempos de cambio de frame (carga y lectura de buffer)
t_lr = interp1(t_er, (1:N)*T_frame, t_frameChange, 'previous', 0); % Tiempo de señal reproducido hasta el momento en cada cambio de buffer

% Tiempo que pasa entre un cambio de de frame y el momento en que ya hay
% suficiente señal en la cola de buffer como para reproducir otro
t_enough = interp1([0; t_lq], [0; t_bufferQueueLoad], t_lr + T_frame, 'next', Inf); % Tiempos en los que ya hay suficiente señal en la cola de buffer como para reproducir otro frame
t_waiting = t_enough - t_frameChange;

% Condición 1: La reproducción de un frame no vacío siempre ha de producirse cuando hay
% suficiente señal en la cola de buffer como para reproducirse. Es decir,
% t_waiting ha de ser negativo siempre que una reproducción se inicie
% 
% Condición 2: La reproducción de un frame vacío siempre ha de producirse cuando todavía
% no hay suficiente señal en la cola de buffer como para reproducirse. Es
% decir, t_waiting ha de ser positivo siempre que una reproducción de frame
% vacío (underrun) ocurra
rbFlag = ismember(t_frameChange, t_br); % Índices de cambios de frame en los que se inicia una reproducción
underrunFlag = ~rbFlag; % Índices de cambios de frame en los que se inicia un underrun (frame vacío)

lowLimitDelay = max(t_waiting(rbFlag)); % Límite inferior de retardo que puede añadirse para que se cumpla la condición de causalidad (condición 1)
upLimitDelay = min(t_waiting(underrunFlag)); % Límite superior de retardo que puede añadirse para que se cumpla la condición 2
assert(upLimitDelay >= lowLimitDelay, 'Error!!!')
delayLimits = [lowLimitDelay, upLimitDelay];
% differ = upLimitDelay - lowLimitDelay; % Margen de indeterminación
% taddDelay = max(0, lowLimitDelay); % Retardo que hay que añadirse al tiempo de inicio de reproducción actual
% tSubstractDelay = max(0, -upLimitDelay); % Retardo que hay que sustraer al tiempo de inicio de reproducción actual

% Representar
% La señal acumulada en la cola de buffer evoluciona discontinuamente. Hay
% que representarlo a escalones.
x = kron(t_bufferQueueLoad, [1; 1]);
base = 0; % Cantidad de señal que había antes de la primera carga de cola
y = reshape([[base; t_lq(1:end-1)], t_lq]', 2*N, 1);

fig = figure;
ax = axes(fig);
plot(ax, t_frameChange, t_lr, x, y);

end

