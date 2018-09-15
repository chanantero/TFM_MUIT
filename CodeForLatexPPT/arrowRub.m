function patchObject = arrowRub( origin, ending, varargin )

p = inputParser;

addParameter(p, 'headLength', []) % Data units
addParameter(p, 'baseAngle', 90) % Degrees
addParameter(p, 'tipAngle', 16) % Degrees
addParameter(p, 'lineWidth', 0)
addParameter(p, 'headNormal', []) % Vector. Normal direction to the head. At least, it should be different from the arrow direction
addParameter(p, 'axes', []);

parse(p, varargin{:})

headLength = p.Results.headLength;
baseAngle = p.Results.baseAngle;
tipAngle = p.Results.tipAngle;
lineWidth = p.Results.lineWidth;
headNormal = p.Results.headNormal;
ax = p.Results.axes;

% Calculate arrow points in a stardardized plane. XY plane and the arrow
% pointing in the positive X direction.
arrowLength = norm(ending - origin); % Euclidean distance: sqrt(sum((ending - origin).^2));
if isempty(headLength)
    headLength = arrowLength * 0.1;
end
tipHeadNormPos = [arrowLength - headLength, tand(tipAngle)*headLength];
baseHeadNormPos = [tipHeadNormPos(1) + tipHeadNormPos(2)/tand(baseAngle), lineWidth/2];

pointsAux1 = [...
    0, lineWidth/2;
    baseHeadNormPos;
    tipHeadNormPos];
    
pointsAux2 = flipud(pointsAux1);
pointsAux2(:, 2) = -pointsAux2(:, 2);

points = [pointsAux1; [arrowLength 0]; pointsAux2];
points = [points, zeros(size(points,1), 1)];

% Transform (rotate) arrow points to the right plane
rotvec = vrrotvec([1 0 0], ending - origin); % [xAxis, yAxis, zAxis, angleOfRotationCounterClockWise]
rotmat = vrrotvec2mat(rotvec);
points = points*rotmat';

% Add rotation in order to orient the head adequately
if ~isempty(headNormal)
headNormalAfterRotation = [0 0 1]*rotmat';
headNormalUser = headNormal - (ending - origin)*dot(ending - origin, headNormal)/arrowLength^2;
headrotvec = vrrotvec(headNormalAfterRotation, headNormalUser); % [xAxis, yAxis, zAxis, angleOfRotationCounterClockWise]
headrotmat = vrrotvec2mat(headrotvec);
points = points*headrotmat';
end

% Translate it to the adequate position
points = points + repmat(origin, [size(points,1), 1]);

if isempty(ax)
ax = gca;
end

patchObject = patch(ax, 'XData', points(:,1), 'YData', points(:,2), 'ZData', points(:,3));
end

