function [ x, y, alfa ] = octogon(d, nb, nd, nl, betabd)

% d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
% nb = 8; % Bottom and upper sides of the octogon
% nd = 8; % Diagonal sides of the octogon (4 sides)
% nl = 24; % Lateral side of the octogon (2 sides)
% betabd = 45; % Deviation angle between bottom/upper and diagonal sides

N = nb*2 + nl*2 + nd*4;
x = zeros(N, 1);
y = zeros(N, 1);
alfa = zeros(N, 1);

v = @(n) linspace(0.5, 0.5 + (n-1), n)*d; % Center of the loudspeakers' front

k = 0;

% Diagonal bottom-left
ind = k+1:k+nd;
x(ind) = flip(v(nd)*cosd(betabd));
y(ind) = v(nd)*sind(betabd);
alfa(ind) = 90 - betabd;
k = k + nd;

% Left lateral
ind = k+1:k+nl;
x(ind) = 0;
y(ind) = v(nl) + nb*d*sind(betabd);
alfa(ind) = 0;
k = k + nl;

% Diagonal top-left
ind = k+1:k+nd;
x(ind) = v(nd)*cosd(betabd);
y(ind) = v(nd)*sind(betabd) + (sind(betabd)*nd*d + nl*d);
alfa(ind) = -90 + betabd;
k = k + nd;

% Top
ind = k+1:k+nb;
x(ind) = v(nb) + cosd(betabd)*nd*d;
y(ind) = 2*sind(betabd)*nd*d + nl*d;
alfa(ind) = -90;
k = k + nb;

% Diagonal top-right
ind = k+1:k+nd;
x(ind) = v(nd)*cosd(betabd) + nb*d + cosd(betabd)*nd*d;
y(ind) = -v(nd)*sind(betabd) + 2*sind(betabd)*nd*d + nl*d;
alfa(ind) = -90 - betabd;
k = k + nd;

% Right lateral
ind = k+1:k+nl;
x(ind) = nb*d + 2*cosd(betabd)*nd*d;
y(ind) = -v(nl) + (sind(betabd)*nd*d + nl*d);
alfa(ind) = 180;
k = k + nl;

% Diagonal bottom-right
ind = k+1:k+nd;
x(ind) = -v(nd)*cosd(betabd) + nb*d + 2*cosd(betabd)*nd*d;
y(ind) = -v(nd)*sind(betabd) + sind(betabd)*nd*d;
alfa(ind) = 180 - betabd;
k = k + nd;

% Bottom
ind = k+1:k+nb;
x(ind) = -v(nb) + nb*d + cosd(betabd)*nd*d;
y(ind) = 0;
alfa(ind) = 90;

end

