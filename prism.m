function [ r, n ] = prism( XLim, YLim, ZLim, XnumPoints, YnumPoints, ZnumPoints )
% n is the outward normal vector. It is not unitary. The magnitude is equal
% to the area of the surface differential.

xVec = (0:XnumPoints-1)'*diff(XLim)/XnumPoints + XLim(1);
yVec = (0:YnumPoints-1)'*diff(YLim)/YnumPoints + YLim(1);
zVec = (0:ZnumPoints-1)'*diff(ZLim)/ZnumPoints + ZLim(1);

dx = xVec(2) - xVec(1);
dy = yVec(2) - yVec(1);
dz = zVec(2) - zVec(1);

% Move to the center of the square
xVec = xVec + dx/2;
yVec = yVec + dy/2;
zVec = zVec + dz/2;

% Inferior face
r1 = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), ones(XnumPoints*YnumPoints, 1)*ZLim(1)];
n1 = repmat([0, 0, -1], size(r1, 1), 1)*(dx*dy);

% Superior face
r2 = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), ones(XnumPoints*YnumPoints, 1)*ZLim(2)];
n2 = repmat([0, 0, 1], size(r2, 1), 1)*(dx*dy);

% Left face
r3 = [XLim(1)*ones(ZnumPoints*YnumPoints, 1), kron(yVec, ones(ZnumPoints, 1)), repmat(zVec, YnumPoints, 1)];
n3 = repmat([-1, 0, 0], size(r3, 1), 1)*(dz*dy);

% Right face
r4 = [XLim(2)*ones(ZnumPoints*YnumPoints, 1), kron(yVec, ones(ZnumPoints, 1)), repmat(zVec, YnumPoints, 1)];
n4 = repmat([1, 0, 0], size(r4, 1), 1)*(dz*dy);

% Rear face
r5 = [kron(xVec, ones(ZnumPoints, 1)), YLim(1)*ones(ZnumPoints*XnumPoints, 1), repmat(zVec, XnumPoints, 1)];
n5 = repmat([0, -1, 0], size(r5, 1), 1)*(dx*dz);

% Forward face
r6 = [kron(xVec, ones(ZnumPoints, 1)), YLim(2)*ones(ZnumPoints*XnumPoints, 1), repmat(zVec, XnumPoints, 1)];
n6 = repmat([0, 1, 0], size(r6, 1), 1)*(dx*dz);

r = [r1; r2; r3; r4; r5; r6];

n = [n1; n2; n3; n4; n5; n6];

end

