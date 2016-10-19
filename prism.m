function [ r, n ] = prism( XLim, YLim, ZLim, XnumPoints, YnumPoints, ZnumPoints )

xVec = linspace(XLim(1), XLim(2), XnumPoints)';
yVec = linspace(YLim(1), YLim(2), YnumPoints)';
zVec = linspace(ZLim(1), ZLim(2), ZnumPoints)';

dx = xVec(2) - xVec(1);
dy = yVec(2) - yVec(1);
dz = zVec(2) - zVec(1);

% Inferior face
r1 = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), ones(XnumPoints*YnumPoints, 1)*ZLim(1)];
n1 = repmat([0, 0, 1], size(r1, 1), 1)*(dx*dy);

% Superior face
r2 = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), ones(XnumPoints*YnumPoints, 1)*ZLim(2)];
n2 = repmat([0, 0, -1], size(r2, 1), 1)*(dx*dy);

% Left face
r3 = [XLim(1)*ones((ZnumPoints-2)*YnumPoints, 1), kron(yVec, ones(ZnumPoints-2, 1)), repmat(zVec(2:end-1), YnumPoints, 1)];
n3 = repmat([1, 0, 0], size(r3, 1), 1)*(dz*dy);

% Right face
r4 = [XLim(2)*ones((ZnumPoints-2)*YnumPoints, 1), kron(yVec, ones(ZnumPoints-2, 1)), repmat(zVec(2:end-1), YnumPoints, 1)];
n4 = repmat([-1, 0, 0], size(r4, 1), 1)*(dz*dy);

% Rear face
r5 = [kron(xVec(2:end-1), ones(ZnumPoints-2, 1)), YLim(1)*ones((ZnumPoints-2)*(XnumPoints-2), 1), repmat(zVec(2:end-1), (XnumPoints-2), 1)];
n5 = repmat([0, 1, 0], size(r5, 1), 1)*(dx*dz);

% Forward face
r6 = [kron(xVec(2:end-1), ones(ZnumPoints-2, 1)), YLim(2)*ones((ZnumPoints-2)*(XnumPoints-2), 1), repmat(zVec(2:end-1), XnumPoints-2, 1)];
n6 = repmat([0, -1, 0], size(r6, 1), 1)*(dx*dz);

r = [r1; r2; r3; r4; r5; r6];

n = [n1; n2; n3; n4; n5; n6];

end

