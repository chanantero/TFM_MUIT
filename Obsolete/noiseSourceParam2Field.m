function field = noiseSourceParam2Field(WFSToolObject, parameters)

parameters = num2cell(parameters);
[x, y, z, amplitude, phase] = deal(parameters{:});

% Set virtual noise source (2nd source) parameters
pos = [x, y, z];
WFSToolObject.amplitude(2) = amplitude;
WFSToolObject.phase(2) = phase;
WFSToolObject.noiseSourcePosition(2, :) = pos;

% WFS calculation
WFSToolObject.WFScalculation(...
    'SourceFilter', 'Loudspeakers',...
    'AcousticPath', 'Theoric',...
    'Grouping', 'AllTogether',...
    'maxAbsoluteValueConstraint', true);

% Simulate
WFSToolObject.simulate();

% Output
field = WFSToolObject.simulField;

end