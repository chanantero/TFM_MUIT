function [r] = square(lim1, lim2, N1, N2, u, theta, desp)

% Input arguments:
% - lim1. Limits along X of the original plane (first quadrant, XY plane)
% - lim2. Limits along Y of the original plane (first quadrant, XY plane)
% - N1. Number of points along X of the original plane
% - N2. Number of points along Y of the original plane
% - desp. Displacement of the original plane. 3 element vector
% - u. Axis of the rotation
% - theta. Rotation on the vector u in clocwise direction
%
% Output arguments:
% - r. Nx3 matrix where each row is a point
% - n. Normal to the surface element

% Create original square
vecX = linspace(lim1(1), lim1(2), N1)';
vecY = linspace(lim2(1), lim2(2), N2)';

rBot = [vecX, lim2(1)*ones(N1, 1)]; % Bottom side
rTop = [vecX, lim2(2)*ones(N1, 1)]; % Top side
rLeft = [lim1(1)*ones(N2-2, 1), vecY(2:end-1)]; % Left side
rRight = [lim1(2)*ones(N2-2, 1), vecY(2:end-1)]; % Right side

rXY = [[rBot; rTop; rLeft; rRight], zeros(N1*2 + N2*2 - 4, 1)];

% Create quaternion for rotation and rotate. The point of rotation is
% [0,0,0]
if ~isempty(u)
    q = [cos(theta/2), u*sin(theta/2)]; % Create quaternion
    r = quatrotate(q, rXY); % Rotate
else
    r = rXY;
end

% Move it
if ~isempty(desp)
    r = r + repmat(desp, size(rXY, 1), 1);
end

end