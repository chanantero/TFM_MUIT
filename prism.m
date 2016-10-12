function [ r ] = prism( XLim, YLim, ZLim, XnumPoints, YnumPoints, ZnumPoints )

xVec = linspace(XLim(1), XLim(2), XnumPoints)';
yVec = linspace(YLim(1), YLim(2), YnumPoints)';
zVec = linspace(ZLim(1), ZLim(2), ZnumPoints)';

% Inferior face
r1 = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), ones(XnumPoints*YnumPoints, 1)*ZLim(1)];

% Superior face
r2 = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), ones(XnumPoints*YnumPoints, 1)*ZLim(2)];

% Left face
r3 = [XLim(1)*ones(ZnumPoints*(YnumPoints-2), 1), kron(yVec, ones(ZnumPoints-2, 1)), repmat(zVec(2:end-1), YnumPoints, 1)];

% Right face
r4 = [XLim(2)*ones(ZnumPoints*(YnumPoints-2), 1), kron(yVec, ones(ZnumPoints-2, 1)), repmat(zVec(2:end-1), YnumPoints, 1)];

% Rear face
r5 = [kron(xVec(2:end-1), ones(ZnumPoints-2, 1)), YLim(1)*ones((ZnumPoints-2)*(XnumPoints-2), 1), repmat(zVec(2:end-1), (XnumPoints-2), 1)];

% Forward face
r6 = [kron(xVec(2:end-1), ones(ZnumPoints-2, 1)), YLim(2)*ones((ZnumPoints-2)*(XnumPoints-2), 1), repmat(zVec(2:end-1), XnumPoints-2, 1)];

r = [r1; r2; r3; r4; r5; r6];


end

