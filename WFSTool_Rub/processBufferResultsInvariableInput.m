function [ fig, delayLimits ] = processBufferResultsInvariableInput( t_bufferQueueLoad, t_underrun, T_frame)
% Invariable input size

N = numel(t_bufferQueueLoad); % N�mero de frames cargados. numel(t_bufferQueueLoad) == numel(t_underrun)

% Buffer theory
% t_bufferQueueLoad: Tiempos de carga de frame en la cola de buffer
t_lq = T_frame*(1:N)'; % Tiempo de se�al acumulado en la cola de buffer en los tiempos de carga en la cola de buffer (despu�s de la carga)

t_delay = 0;
t_br = t_delay + cumsum([t_underrun(1); t_underrun(2:end) + T_frame]); % Tiempos de inicio de reproducci�n
t_er = t_br + T_frame; % Tiempos de fin de reproducci�n

t_frameChange = (t_br(1):T_frame:t_er(end))'; % Tiempos de cambio de frame (carga y lectura de buffer)
t_lr = interp1(t_er, (1:N)*T_frame, t_frameChange, 'previous', 0); % Tiempo de se�al reproducido hasta el momento en cada cambio de buffer

% Tiempo que pasa entre un cambio de de frame y el momento en que ya hay
% suficiente se�al en la cola de buffer como para reproducir otro
t_enough = interp1([0; t_lq], [0; t_bufferQueueLoad], t_lr + T_frame, 'next', Inf); % Tiempos en los que ya hay suficiente se�al en la cola de buffer como para reproducir otro frame
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
delayLimits = [lowLimitDelay, upLimitDelay];
% differ = upLimitDelay - lowLimitDelay; % Margen de indeterminaci�n
% taddDelay = max(0, lowLimitDelay); % Retardo que hay que a�adirse al tiempo de inicio de reproducci�n actual
% tSubstractDelay = max(0, -upLimitDelay); % Retardo que hay que sustraer al tiempo de inicio de reproducci�n actual

% Representar
% La se�al acumulada en la cola de buffer evoluciona discontinuamente. Hay
% que representarlo a escalones.
x = kron(t_bufferQueueLoad, [1; 1]);
base = 0; % Cantidad de se�al que hab�a antes de la primera carga de cola
y = reshape([[base; t_lq(1:end-1)], t_lq]', 2*N, 1);

fig = figure;
ax = axes(fig);
plot(ax, t_frameChange, t_lr, x, y);

end

