function [ x, y ] = octogon(d, nb, nd, nl, betabd)

% d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
% nb = 8; % Bottom and upper sides of the octogon
% nd = 8; % Diagonal sides of the octogon (4 sides)
% nl = 24; % Lateral side of the octogon (2 sides)
% betabd = 45; % Deviation angle between bottom/upper and diagonal sides

N = nb*2 + nl*2 + nd*4;
x = zeros(N, 1);
y = zeros(N, 1);

v = @(n) linspace(0.5, 0.5 + (n-1), n)*d; % Center of the loudspeakers' front

k = 0;

x(k+1:k+nd) = flip(v(nd)*cosd(betabd));
y(k+1:k+nd) = v(nd)*sind(betabd);
k = k + nd;

x(k+1:k+nl) = 0;
y(k+1:k+nl) = v(nl) + nb*d*sind(betabd);
k = k + nl;

x(k+1:k+nd) = v(nd)*cosd(betabd);
y(k+1:k+nd) = v(nd)*sind(betabd) + (sind(betabd)*nd*d + nl*d);
k = k + nd;

x(k+1:k+nb) = v(nb) + cosd(betabd)*nd*d;
y(k+1:k+nb) = 2*sind(betabd)*nd*d + nl*d;
k = k + nb;

x(k+1:k+nd) = v(nd)*cosd(betabd) + nb*d + cosd(betabd)*nd*d;
y(k+1:k+nd) = -v(nd)*sind(betabd) + 2*sind(betabd)*nd*d + nl*d;
k = k + nd;

x(k+1:k+nl) = nb*d + 2*cosd(betabd)*nd*d;
y(k+1:k+nl) = -v(nl) + (sind(betabd)*nd*d + nl*d);
k = k + nl;

x(k+1:k+nd) = -v(nd)*cosd(betabd) + nb*d + 2*cosd(betabd)*nd*d;
y(k+1:k+nd) = -v(nd)*sind(betabd) + sind(betabd)*nd*d;
k = k + nd;

x(k+1:k+nb) = -v(nb) + nb*d + cosd(betabd)*nd*d;
y(k+1:k+nb) = 0;

end

