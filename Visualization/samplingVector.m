function vector = samplingVector( positions, minimum, maximum, maximumSpacing, minNumIntervBetween2Pos)

positions = sort(positions);
numSegments = numel(positions) + 1;
segments = cell(1, numSegments);

for s = 1:numSegments
    if(s < numSegments)
        highValue = positions(s);
    else
        highValue = maximum;
    end
    if(s > 1)
        lowValue = positions(s - 1);
    else
        lowValue = minimum;
    end
    
    step = min((highValue - lowValue)/ceil((highValue - lowValue)/maximumSpacing), (highValue - lowValue)/minNumIntervBetween2Pos);
    segment = lowValue:step:highValue;
    if(step==0) % isempty(segment)
        segment = lowValue;
    end
    
    segments{s} = segment;
end

vector = unique(cell2mat(segments));

end

