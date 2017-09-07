function acceptable = isAccepted( realDoAs, calcDoAs, indRealDoAs, indCalcDoAs, maxDev )

realDoAs = realDoAs(:);
calcDoAs = calcDoAs(:);
indRealDoAs = indRealDoAs(:);
indCalcDoAs = indCalcDoAs(:);

if(nargin < 5)
    maxDev = 10;
end

acceptable = false;
if(numel(calcDoAs) == numel(realDoAs))
    if(abs(realDoAs(indRealDoAs) - calcDoAs(indCalcDoAs)) <= maxDev)
        acceptable = true;
    end
end

end

