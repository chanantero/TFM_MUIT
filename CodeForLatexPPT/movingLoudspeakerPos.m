function [ pos ] = movingLoudspeakerPos( shape, varargin )

switch shape
    case 'eight'
        p = inputParser;
        
        addParameter(p, 'numPoints', 100);
        addParameter(p, 'radius', 0.8);
        addParameter(p, 'refPos', [0 0]);
        
        parse(p, varargin{:})
        
        numPoints = p.Results.numPoints;
        radius = p.Results.radius;
        refPos = p.Results.refPos;
        
        refPosX = refPos(1);
        refPosY = refPos(2);
        
        % Virtual Position
        numPosCircle = floor(numPoints/2);
        alpha = (0:numPosCircle-1)/numPosCircle * 2*pi;
        
        centreX1 = refPosX;
        centreY1 = refPosY + radius;
        circleX1 = centreX1 + radius*cos(alpha - pi/2);
        circleY1 = centreY1 + radius*sin(alpha - pi/2);
        
        centreX2 = refPosX;
        centreY2 = refPosY - radius;
        circleX2 = centreX2 + radius*cos(-alpha + pi/2);
        circleY2 = centreY2 + radius*sin(-alpha + pi/2);
        
        pos = [circleX1' circleY1'; circleX2' circleY2'];
        pos = [pos, zeros(size(pos, 1), 1)];
        
    case 'spiral'
        p = inputParser;
        
        addParameter(p, 'numPoints', 100);
        addParameter(p, 'numRounds', 1);
        addParameter(p, 'radius', 0.8);
        addParameter(p, 'refPos', [0 0]);
        
        parse(p, varargin{:})
        
        numPoints = p.Results.numPoints;
        numRounds = p.Results.numRounds;
        radius = p.Results.radius;
        refPos = p.Results.refPos;
        
        refPosX = refPos(1);
        refPosY = refPos(2);
        
        rStep = 2/numPoints;
        r = (1:-rStep:-(1-rStep))';
        r(r<0) = -r(r<0);
        r = r*radius;
        theta = linspace(0, numRounds*2*pi*2, numPoints)';
        
        [x, y] = pol2cart(theta, r);
        
        x = x + refPosX;
        y = y + refPosY;
        
        pos = [x, y];
        
        pos = [pos, zeros(size(pos, 1), 1)];
        
end


end

