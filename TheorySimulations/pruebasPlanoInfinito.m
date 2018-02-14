%% Plano infinito

% Plano infinito. Cada diferencial de superficie es una fuente monopolo
% En un punto situado a una distancia z, del plano, sumamos la contribución
% de todos los puntos. Esto lo calcularé con una integral, no con una
% simulación de muchos monopolos
k = 1;
field_r = @(z, r) exp(-1i*k*sqrt(z^2 + r.^2))./sqrt(z^2 + r.^2)*2*pi.*r;

r = 1e4:0.01:1e4+10;
ax = axes(figure, 'NextPlot', 'Add');
plot(ax, r, real(field_r(1, r)));
plot(ax, r, imag(field_r(1, r)));

z = 0.1:0.01:50;
N = length(z);
field = zeros(N, 1);
for k = 1:N
    fAux = @(r) field_r(z(k), r);
    field(k) = integral(fAux, 0, 1e4);
end
plot(z, field)


% Método más bruto. Fuentes en el plano x y, punto de recepción situado
% en el eje z.
field_xy = @(x, y, z) exp(-1i*k*sqrt(x.^2 + y.^2 + z.^2))./sqrt(x.^2 + y.^2 + z.^2);

z = 1:0.1:50;
L = 40; step = 0.01;
xLim = [-L/2, L/2];
yLim = [-L/2, L/2];
dx = step;
dy = step;
x = xLim(1):dx:xLim(2);
y = yLim(1):dy:yLim(2);
[X, Y] = ndgrid(x, y);

N = length(z);
field = zeros(N, 1);
n = 0;
for k = 1:N
    fprintf(repmat('\b', 1, n));
    msg = sprintf('%d/%d', k, N);
    fprintf(msg);
    n = numel(msg);
    
%     fAux = @(x, y) field_xy(x, y, z(k));
%     field(k) = integral2(fAux, xLim(1), xLim(2), yLim(1), yLim(2));
    
    field(k) = sum(field_xy(X(:), Y(:), z(k)))*dx*dy;
end

ax = axes(figure);
plot(z, real(field), z, imag(field))

%% Concatenación de hilos infinitos
% Podemos considerar a un plano infinito como la concatenación de hilos
% infinitos paralelos

% Contribución del hilo paralelo al eje X, situado en el plano z = 0, en la
% posición y, al punto situado en el eje Z ([0 0 z]).
k = 1;
z = 0.01:0.05:50;
f = @(y, z) -pi*1i*besselh(0, 2, k*sqrt(z.^2 + y.^2));

N = length(z);
field = zeros(N, 1);
n = 0;
for k = 1:N
    fprintf(repmat('\b', 1, n));
    msg = sprintf('%d/%d', k, N);
    fprintf(msg);
    n = numel(msg);
    
    fAux = @(y) f(y, z(k));
    field(k) = integral(fAux, -1000, 1000);
end

ax = axes(figure);
plot(ax, z, real(field), z, imag(field))